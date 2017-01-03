
open Std

let suffix = "(.txt)?$"

(* The nuc/gen .txt are all optional at the moment, so we need the '$'
   to bind to the end! *)
let filename_regex =
  Re_posix.compile_pat ("([^/]+)/([^/]+)_[^/]+(_(nuc|gen))?" ^ suffix)

let allele_to_hla_class s =
  if String.get s 0 = 'D' then II else I

let parse (fname, re_group) =
  let sample = Re.Group.get re_group 2 in
  let sample =
    if String.contains sample '_' then
      String.sub sample 0 (String.index sample '_')
    else
      sample
  in
  let ic = open_in fname in
  let rec loop acc =
    try
      let line = input_line ic in
      let confidence, allele = Scanf.sscanf line "%s\t%s" (fun c a -> (c, a)) in
      let info =
        let open Info in
        { hla_class = allele_to_hla_class allele
        ; qualifier = ""
        ; allele
        ; confidence = float_of_string_nanable confidence
        ; typer_spec = ""
        }
      in
      loop (info :: acc)
    with End_of_file -> acc
  in
  let ilst =
    try loop []
    with e ->
      Printf.eprintf "dropping %s because of %s\n" fname (Printexc.to_string e);
      []
  in
  sample, ilst

let scan_directory (dir : string) : ((sample * Info.t list) list) =
  let rec loop acc = function
    | [] -> acc
    | dir :: t ->
      let matches, dirs =
        Sys.readdir dir
        |> Array.fold_left ~init:([],[]) ~f:(fun (ms, ds) f ->
              let df = Filename.concat dir f in
              if Sys.is_directory df then
                ms, df :: ds
              else
                match Re.exec_opt filename_regex df with
                | None    -> ms, ds
                | Some m  -> (df, m) :: ms, ds)
      in
      loop (matches @ acc) (dirs @ t)
  in
  loop [] [dir]
  |> List.map ~f:parse
  |> join_by_fst
  |> List.sort ~cmp:compare
