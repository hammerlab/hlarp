
open Std

let suffix = ".HLAgenotype4digits"
let filename_regex = Re_posix.compile_pat ("(.+)-Class(I{1,2})" ^ suffix)

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
          Scanf.sscanf (input_line ic) "%s\t%s\t%s\t%s\t%s"
            (fun _locus allele1 conf1 allele2 conf2 ->
                { hla_class  = cls
                ; allele     = allele1
                ; qualifier  = ""
                ; confidence = float_of_string_nanable conf1
                },
                { hla_class  = cls
                ; allele     = allele2
                ; qualifier  = ""
                ; confidence = float_of_string_nanable conf2
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
  |> join_by_fst
  |> List.sort ~cmp:compare  (* sort by keys aka runs *)

