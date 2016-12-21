(** A tool to normalize HLA typing output.

class,allele,qualifier,confidence,typer_spec,sample
1,HLA-A*02:01,G,0.9,,Patient_1
*)

open Hlarp
open Std

let project_name = "hlarp"

let multiple seqlst optlst athlst prolst do_not_prefix =
  let p typer_name =
    if do_not_prefix then
      fun x -> x
    else
      List.map ~f:(fun (p,l) -> (typer_name ^ p, l))
  in
  [ List.map ~f:Seq2HLA.scan_directory seqlst |> List.map ~f:(p "seq2HLA_")
  ; List.map ~f:OptiType.scan_directory optlst |> List.map ~f:(p "OptiType_")
  ; List.map ~f:Athlates.scan_directory athlst |> List.map ~f:(p "ATHLATES_")
  ; List.map ~f:Prohlatype.scan_directory prolst |> List.map ~f:(p "Prohlatype")
  ]
  |> List.concat
  |> List.concat
  |> Output.out_channel stdout

let compare resolution classes loci max_allele_rows_to_print metrics summary_by
  (* Input *)
  pos_args =
  let classes = match classes with | [] -> None | l -> Some l in
  let loci = match loci with | [] -> None | l -> Some l in
  let nmp  = Compare.nested_maps pos_args in
  let azo  = Compare.analyze_samples ?loci ?resolution ?classes ~metrics nmp in
  Compare.output ?max_allele_rows_to_print ~summary_by stdout azo

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
    let f dir = Seq2HLA.scan_directory dir |> Output.out_channel stdout in
    Term.(const f $ (directory_arg ~tool ~suffix:Seq2HLA.suffix)
        , info tool ~doc:"Parse seq2HLA output")
  in
  let optitype =
    let tool = "optitype" in
    let f dir = OptiType.scan_directory dir |> Output.out_channel stdout in
    Term.(const f $ (directory_arg ~tool ~suffix:OptiType.suffix)
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
    let f dir equal_pairs =
      Athlates.scan_directory ~equal_pairs dir |> Output.out_channel stdout
    in
    Term.(const f
          $ directory_arg ~tool ~suffix:Athlates.suffix
          $ equal_pairs_flag
        , info tool ~doc:"Parse ATHLATES output")
  in
  let prohlatype =
    let tool = "prohlatype" in
    let f dir = Prohlatype.scan_directory dir |> Output.out_channel stdout in
    Term.(const f $ directory_arg ~tool ~suffix:Prohlatype.suffix
        , info tool ~doc:"Parse Prohlatype output")
  in
  let to_directory_arg name arg_lst =
    Arg.(value & opt_all dir []
        & info arg_lst ~docv:"DIR"
            ~doc:(sprintf "Directory to search for %s output (repeatable)." name))
  in
  let seq_arg = to_directory_arg "Seq2HLA" [ "seq2HLA" ] in
  let opt_arg = to_directory_arg "Optitype" [ "optitype" ] in
  let ath_arg = to_directory_arg "ATHLATES" [ "athlates" ] in
  let pro_arg = to_directory_arg "Prohlatype" [ "prohlatype" ] in
  let typers_as_pos_args =
    let is_file f o arg =
      if Sys.file_exists f then
        `Ok o
      else
        `Error (sprintf "not a file %s in %s argument" f arg)
    in
    let is_dir s o arg =
      if Sys.file_exists s && Sys.is_directory s then
        `Ok o
      else
        `Error (sprintf "not a directory %s in %s argument" s arg)
    in
    let tparser s =
      try
        let i = String.index s '=' in
        let n = String.length s in
        let l = String.sub s 0 i in
        let r = String.sub s (i + 1) (n - i - 1) in
        match l with
        | "hlarp-file" -> is_file r (`HlarpFile r) l
        | "seq2HLA"    -> is_dir r (`seq2HLA r) l
        | "optitype"   -> is_dir r (`OptiType r) l
        | "athlates"   -> is_dir r (`ATHLATES r) l
        | "prohlatype" -> is_dir r (`Prohlatype r) l
        | s            -> `Error ("unrecognized typer spec: " ^ s)
      with Not_found ->
        `Error (sprintf "Failed to find '=' in typer specific argument: %s" s)
    in
    let tprinter fmt = function
      | `HlarpFile f  -> Format.fprintf fmt "Hlarp file %s" f
      | `seq2HLA d    -> Format.fprintf fmt "seq2HLA %s" d
      | `OptiType d   -> Format.fprintf fmt "OptiType %s" d
      | `ATHLATES d   -> Format.fprintf fmt "ATHLATES %s" d
      | `Prohlatype d -> Format.fprintf fmt "Prohlatype %s" d
    in
    Arg.(non_empty
        & pos_all (tparser, tprinter) []
        & info
          ~docv:"[seq2HLA|optitype|athlates|prohlatype|hlarp-file]=[DIRECTORY|FILE]"
          ~doc:"Directories or files (as in the last case) prefixed by the \
                name of the typer and an equal sign, to specify where to look \
                for aggregation or comparison HLA data."
          [])
  in
  let multiple =
    let pre_flg =
      Arg.(value & flag & info ["prefix"]
            ~doc:"Do NOT prefix the run information with typer name.")
    in
    Term.(const multiple $ seq_arg $ opt_arg $ ath_arg $ pro_arg $ pre_flg
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
    let loci_arg =
      Arg.(value & opt_all string []
            & info ["l"; "loci"]
                ~doc:"MHC loci along which to partition the similarity analysis: a string prefix.\
                      Specify a prefix that is used to group loci (ex. \"A\", \"B\", \"DRB1\", etc.).\
                      Specify multiple loci to get separate analyses.\
                      This argument will supersede any class grouping arguments.")
    in
    let max_allele_rows_to_print_arg =
      Arg.(value & opt (some int) None
          & info ["max-allele-rows-to-print"]
              ~docv:"NON NEGATIVE INTEGER"
              ~doc:"The number of per sample per comparison alleles to \
                    display. Defaults to 0")
    in
    let metrics_flag =
      Arg.(value & vflag_all [`Jaccard]
        [ `Jaccard, info ~doc:"Jaccard similarity" ["jaccard"]
        ; `KLDiv,   info ~doc:"Kullbeck Leibler divergence" ["kldiv"]
        ])
    in
    let summary_by_arg =
      Arg.(value & vflag `Mean
        [ `Mean,    info ~doc:"Report the mean across samples by comparison." ["mean-summary"]
        ; `Median,  info ~doc:"Report the median across samples by comparison." ["median-summary"]
        ])
    in
    Term.(const compare $ resolution_arg $ classes_arg $ loci_arg
            $ max_allele_rows_to_print_arg
            $ metrics_flag
            $ summary_by_arg
            (* Input *)
            $ typers_as_pos_args
        , info "compare"
            ~doc:"Scan multiple directories (of possibly different formats) and compare the results after aggregating on a per run basis.")
  in
  let cmds = [seq2HLA; optitype; athlates; prohlatype; multiple; compare] in
  match Term.eval_choice help_cmd cmds with
  | `Ok () -> ()
  | `Error _ -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
