(** A tool to normalize HLA typing output.

class,allele,qualifier,confidence
1,HLA-A*02:01,G,0.9
*)

open Nonstd
open Printf

open Hlarp

let project_name = "hlarp"

let multiple seqlst optlst athlst do_not_prefix =
  let p typer_name =
    if do_not_prefix then
      fun x -> x
    else
      List.map ~f:(fun (p,l) -> (typer_name ^ p, l))
  in
  [ List.map ~f:Seq2HLA.scan_directory seqlst |> List.map ~f:(p "seq2HLA_")
  ; List.map ~f:OptiType.scan_directory optlst |> List.map ~f:(p "OptiType_")
  ; List.map ~f:Athlates.scan_directory athlst |> List.map ~f:(p "ATHLATES_")
  ]
  |> List.concat
  |> List.concat
  |> Output.out_channel stdout

let compare resolution classes seqlst optlst athlst =
  let classes = match classes with | [] -> None | l -> Some l in
  Compare.nested_maps seqlst optlst athlst
  |> Compare.output ?resolution ?classes stdout

let () =
  let open Cmdliner in
  let version = "0.0.0" in
  let help_cmd =
    let doc = "HLA typing output normalizer." in
    let bug =
      sprintf "Browse and report new issues at <https://github.com/hammerlab/%s>."
        project_name
    in
    let man =
      [ `S "AUTHORS"
      ; `P "Leonid Rozenberg <leonidr@gmail.com>"
      ; `Noblank
      ; `S "BUGS"
      ; `P bug
      ]
    in
    Term.(ret (pure (`Help (`Plain, None)))
        , info project_name ~version ~doc ~man)
  in
  let directory_arg ~tool ~suffix =
    let docv = "directory" in
    let doc = sprintf
        "Directory to find %s output. %s will search for files with \"%s\" \
         filename suffix. Defaults to looking in the current dir."
          tool project_name suffix
    in
    Arg.(required
         & pos 0 (some dir) (Some ".")
         & info ~doc ~docv [])
  in
  let seq2HLA =
    let tool = "seq2HLA" in
    Term.(const (fun dir ->
            Seq2HLA.scan_directory dir
            |> Output.out_channel stdout)
          $ (directory_arg ~tool ~suffix:Seq2HLA.suffix)
        , info tool ~doc:"Parse seq2HLA output")
  in
  let optitype =
    let tool = "optitype" in
    Term.(const (fun dir ->
            OptiType.scan_directory dir
            |> Output.out_channel stdout)
          $ (directory_arg ~tool ~suffix:OptiType.suffix)
        , info tool ~doc:"Parse OptiType output")
  in
  let athlates =
    let tool = "athlates" in
    let equal_pairs_flag =
      Arg.(value
           & vflag `MostLikelyPair
              [ `MostLikelyPair, info ["most_likely_pair"]
                  ~doc:"For each typing file report the most frequent two alleles (default)."
              ; `ReportAll,      info ["all_pairs"]
                  ~doc:"For each typing file report all pairs (will create duplicates)."
              ; `FirstPair,      info ["first_pair"]
                  ~doc:"For each typing file report only the first pair (may create duplicates)."
              ; `Unique,         info ["unique"]
                  ~doc:"For each typing file report only the unique alleles."
              ]
      )
    in
    Term.(const (fun dir equal_pairs ->
            Athlates.scan_directory ~equal_pairs dir
            |> Output.out_channel stdout)
          $ directory_arg ~tool ~suffix:Athlates.suffix
          $ equal_pairs_flag
        , info tool ~doc:"Parse ATHLATES output")
  in
  let to_directory_arg name arg_lst =
    Arg.(value & opt_all dir []
        & info arg_lst ~docv:"DIR"
            ~doc:(sprintf "Directory to search for %s output (repeatable)." name))
  in
  let seq_arg = to_directory_arg "Seq2HLA" ["s"; "seq2HLA"] in
  let opt_arg = to_directory_arg "Optitype" ["o"; "optitype"] in
  let ath_arg = to_directory_arg "ATHLATES" ["a"; "athlates"] in
  let multiple =
    let pre_flg =
      Arg.(value & flag & info ["prefix"]
            ~doc:"Do NOT prefix the run information with typer name.")
    in
    Term.(const multiple $ seq_arg $ opt_arg $ ath_arg $ pre_flg
        , info "multiple"
            ~doc:"Scan multiple directories (of possilby different formats) and report aggregate results.")
  in
  let compare =
    let resolution_arg =
      let one_to_four_converter =
        let parser_ = function
          | "1" | "2" | "3" | "4" as r -> `Ok (int_of_string r)
          | s -> `Error (sprintf "Not in [1,4]: %s" s)
        in
        let printer f = Format.fprintf f "%d" in
        parser_, printer
      in
      Arg.(value & opt (some one_to_four_converter) None
            & info ["r"; "resolution"]
              ~doc:"MHC allele resolution to use for the analysis. Must be an integer between 1 and 4.")
    in
    let classes_arg =
      let one_or_two_converter =
        let parser_ = function
          | "1" -> `Ok I | "2" -> `Ok II
          | s -> `Error (sprintf "Not in [1,2]: %s" s)
        in
        let printer f c = Format.fprintf f "%s" (hla_class_to_string c) in
        parser_, printer
      in
      Arg.(value & opt_all one_or_two_converter []
            & info ["c"; "class"]
                ~doc:"MHC class along which to partition the similarity analysis: Must be 1 or 2.\
                        Specify multiple classes to get separate analysis.")
    in
    Term.(const compare $ resolution_arg $ classes_arg $ seq_arg $ opt_arg $ ath_arg
        , info "compare"
            ~doc:"Scan multiple directories (of possibly different formats) and compare the results after aggregating on a per run basis.")
  in
  let cmds = [seq2HLA; optitype; athlates; multiple; compare] in
  match Term.eval_choice help_cmd cmds with
  | `Ok () -> ()
  | `Error _ -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
