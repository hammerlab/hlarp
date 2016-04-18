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

  let parse fname =
    let cls_match = Re.exec filename_regex (Filename.basename fname) in
    let run = Re.Group.get cls_match 1 in
    let cls = to_hla_class (Re.Group.get cls_match 2) in
    let ic = open_in fname in
    let hdr = input_line ic in
    if hdr <> "#Locus\tAllele 1\tConfidence\tAllele 2\tConfidence" then
      begin
        close_in ic;
        failwith ("Unsupported header row: " ^ hdr)
      end
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

  let scan_directory dir =
    Sys.readdir dir
    |> Array.to_list
    |> List.map ~f:(Re.matches glob_regex)
    |> List.concat
    |> List.dedup
    |> List.map ~f:(fun f -> parse (Filename.concat dir f))
    |> List.fold_left ~f:(fun acc (run, inf) ->
        match StringMap.find run acc with
        | exception Not_found -> StringMap.add run inf acc
        | clst -> StringMap.add run (clst @ inf) acc)
        ~init:StringMap.empty
    |> StringMap.bindings
    |> List.sort ~cmp:compare  (* sort by keys aka runs *)

end (* Seq2HLA *)

module OptiType = struct

  let suffix = "_result.tsv"

  (* The file is 2 directories down:
      1st is the 'run' name.
      2nd is the date, since there can be multiple datetimes we'll need some
        logic to combine them.
    We'll hard-code Unix forward slashes. *)
  let filename_regex = Re_posix.compile_pat ("([^/]+)/([^/]+)/[^/]+" ^ suffix)

  let parse fname =
    let grp = Re.all filename_regex fname |> List.hd_exn in
    let ri1 = Re.Group.get grp 1 in
    let ri2 = Re.Group.get grp 2 in
    let ic  = open_in fname in
    let hdr = input_line ic in
    if hdr <> "\tA1\tA2\tB1\tB2\tC1\tC2\tReads\tObjective" then
      begin
        close_in ic;
        failwith ("Unsupported header row: " ^ hdr)
      end
    else
      let info allele =
        { hla_class = I ; allele ; qualifier = ""; confidence = nan}
      in
      let lst =
        Scanf.sscanf (input_line ic)
          "0\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f"
          (fun a1 a2 b1 b2 c1 c2 _reads _objective ->
             [ info a1 ; info a2 ; info b1 ; info b2 ; info c1 ; info c2 ])
      in
      close_in ic;
      (ri1, ri2), lst

  let glob_regex = Re.compile (Re_glob.glob ("*" ^ suffix))

  let scan_directory ?(combine_multiday=`TakeLast) dir =
    let rec loop acc = function
      | [] -> acc
      | dir :: t ->
        let dirs, matches =
          Sys.readdir dir         (* Thankfully doesn't return current dir, '.' *)
          |> Array.to_list
          |> List.fold_left ~init:([],[])
                ~f:(fun (d,m) f->
                    let df = Filename.concat dir f in
                    if Sys.is_directory df then
                      df :: d, m
                    else if Re.execp glob_regex df then
                      d, df::m
                    else
                      d, m)
        in
        loop (matches @ acc) (dirs @ t)
    in
    loop [] [dir]
    |> List.map ~f:parse
    |> List.sort ~cmp:compare (* Sort by keys aka-runs *)
    |> (function              (* Do some deduping of the multi-date things *)
        | [] -> []
        | (((r1,_r2),_) :: _) as lst ->
          match combine_multiday with
          | `TakeLast ->
              List.fold_left ~f:(fun (prev_r1, prev_res, acc) ((r1, _r2), lst) ->
                if prev_r1 = r1 then
                  (r1, lst, acc)
                else
                  (r1, lst, (prev_r1, prev_res) :: acc))
              ~init:("not r1" ^ r1, [], [])
              lst
            |> (fun (_, _, acc) -> acc))

end (* OptiType *)

module Output = struct

  let hla_class_to_string = function | I -> "1" | II -> "2"

  let nan_is_empty f = if f <> f then "" else sprintf "%f" f

  let out_channel oc run_assoc =
    fprintf oc "class,allele,qualifier,confidence,run\n";
    List.iter run_assoc ~f:(fun (run, lst) ->
      List.sort ~cmp:compare lst
      |> List.iter ~f:(fun {hla_class; allele; qualifier; confidence} ->
        fprintf oc "%s,%s,%s,%s,%s\n"
          (hla_class_to_string hla_class)
          allele qualifier (nan_is_empty confidence) run))

end (* Output *)

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
  let directory_arg ~tool ~suffix =
    let doc =
      sprintf "Directory to find %s output.\
                      \n%s will search for files with %s suffix.\
                      Defaults to looking in the current dir."
              tool project_name suffix
    in
    Arg.(value & pos 0 dir "." & info ~doc [])
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
  let cmds = [seq2HLA; optitype] in
  match Term.eval_choice help_cmd cmds with
  | `Ok () -> ()
  | `Error _ -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
