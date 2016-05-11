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

let compare resolution seqlst optlst athlst =
  Compare.nested_maps seqlst optlst athlst
  |> Compare.output ?resolution stdout

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
  let seq_arg = Arg.(value & opt_all dir [] & info ~doc:"Seq2HLA dir" ["s"; "seq2HLA"]) in
  let opt_arg = Arg.(value & opt_all dir [] & info ~doc:"Optitype dir" ["o"; "optitype"]) in
  let ath_arg = Arg.(value & opt_all dir [] & info ~doc:"ATHLATES dir" ["a"; "athlates"]) in
  let multiple =
    let pre_flg = Arg.(value & flag & info ~doc:"Do NOT prefix the run information with typer name." ["prefix"]) in
    Term.(const multiple $ seq_arg $ opt_arg $ ath_arg $ pre_flg
        , info "multiple" ~doc:"Multiple")
  in
  let compare =
    let resolution_flag =
      let bounded_int_converter =
        let parser_ = function
          | "1" | "2" | "3" | "4" as r -> `Ok (int_of_string r)
          | s -> `Error (sprintf "Not in [1,4]: %s" s)
        in
        let printer f = Format.fprintf f "%d" in
        parser_, printer
      in
      Arg.(value & opt (some bounded_int_converter) None
            & info ["r"; "resolution"]
              ~doc:"MHC allele resolution to use for the analysis. Must be an integer between 1 and 4.")
    in
    Term.(const compare $ resolution_flag $ seq_arg $ opt_arg $ ath_arg
        , info "compare" ~doc:"Compare")
  in
  let cmds = [seq2HLA; optitype; athlates; multiple; compare] in
  match Term.eval_choice help_cmd cmds with
  | `Ok () -> ()
  | `Error _ -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
