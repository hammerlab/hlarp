(** A tool to normalize HLA typing output. 

class,allele,qualifier,confidence
1,HLA-A*02:01,G,0.9
*)


open Nonstd
open Printf
module StringMap = Map.Make (struct type t = string let compare = compare end) 

type hla_class = | I | II

type info =
  { hla_class   : hla_class
  ; allele      : string
  ; qualifier   : string
  ; confidence  : float
  }

module Seq2HLA = struct

  let suffix = "HLAgenotype4digits"
  let filename_regex = Re_posix.compile_pat ("(.+)-Class(I{1,2})." ^ suffix)

  let to_hla_class = function
    | "I" -> I
    | "II" -> II
    | s -> raise (invalid_arg ("Unrecognized HLA class : " ^ s))

  let parse_seq2HLA_file fname =
    let cls_match = Re.exec filename_regex (Filename.basename fname) in
    let run = Re.Group.get cls_match 1 in
    let cls = to_hla_class (Re.Group.get cls_match 2) in
    let ic = open_in fname in
    let hdr = input_line ic in
    if hdr <> "#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence" then
      failwith ("Unsupported header row: " ^ hdr)
    else
      let rec loop acc =
        try
          let a1, a2 = 
            Scanf.sscanf (input_line ic) "%s\t%s\t%f\t%s\t%f"
              (fun _locus allele1 conf1 allele2 conf2 ->
                 { hla_class  = cls
                 ; allele     = allele1
                 ; qualifier  = ""
                 ; confidence = conf1
                 },
                 { hla_class  = cls
                 ; allele     = allele2
                 ; qualifier  = ""
                 ; confidence = conf2
                 })
          in
          loop (a2 :: a1 :: acc)
        with End_of_file ->
          close_in ic;
          (run, acc)
      in
      loop []

  let glob_regex = Re.compile (Re_glob.glob ("*" ^ suffix))

  let parse_seq2HLA_dir dir =
    Sys.readdir dir
    |> Array.to_list
    |> List.map ~f:(Re.matches glob_regex)
    |> List.concat
    |> List.dedup
    |> List.map ~f:(fun f -> parse_seq2HLA_file (Filename.concat dir f))
    |> List.fold_left ~f:(fun acc (run, inf) ->
        match StringMap.find run acc with
        | exception Not_found -> StringMap.add run inf acc
        | clst -> StringMap.add run (clst @ inf) acc)
        ~init:StringMap.empty
    |> StringMap.bindings
    |> List.sort ~cmp:compare  (* sort by keys aka runs *)

end

module Output = struct

  let hla_class_to_string = function | I -> "1" | II -> "2"

  let out_channel oc run_assoc =
    fprintf oc "class,allele,qualifier,confidence,run\n";
    List.iter run_assoc ~f:(fun (run, lst) ->
      List.sort ~cmp:compare lst 
      |> List.iter ~f:(fun {hla_class; allele; qualifier; confidence} ->
        fprintf oc "%s,%s,%s,%f,%s\n"
          (hla_class_to_string hla_class)
          allele qualifier confidence run))

end

let project_name = "hlarp"

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
  let directory_arg =
    let doc =
      sprintf "Directory to find seq2HLA output.\
                      \n%s will search for files with %s suffix.\
                      Defaults to looking in the current dir."
                project_name Seq2HLA.suffix
    in
    Arg.(value & pos 0 dir "." & info ~doc [])
  in
  let seq2HLA =
    Term.(const (fun dir ->
            Seq2HLA.parse_seq2HLA_dir dir
            |> Output.out_channel stdout)
          $ directory_arg
        , info "seq2HLA" ~doc:"Parse seq2HLA output")
  in
  let cmds = [seq2HLA] in
  match Term.eval_choice help_cmd cmds with
  | `Ok () -> ()
  | `Error _ -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
