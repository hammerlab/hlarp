open Nonstd
open Printf

module StringMap = Map.Make (struct type t = string let compare = compare end)

let join_by_fst lst =
  List.fold_left lst ~f:(fun acc (run, inf) ->
    match StringMap.find run acc with
    | exception Not_found -> StringMap.add run inf acc
    | clst -> StringMap.add run (clst @ inf) acc)
        ~init:StringMap.empty
  |> StringMap.bindings


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
    |> join_by_fst
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
          |> Array.fold_left ~init:([],[]) ~f:(fun (d,m) f->
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

module Athlates = struct

  let suffix = ".typing.txt"

  let filename_regex = Re_posix.compile_pat ("([^/]+)/[^/]+" ^ suffix)

  let header_start = "------------ Inferred Allelic Pairs -------------"

  let allele_to_hla_class s =
    if String.get s 0 = 'D' then II else I

  let parse ?(equal_pairs=`ReportAll) (fname, re_group) =
    let run = Re.Group.get re_group 1 in
    let ic = open_in fname in
    let rec loop found_header acc =
      try
        let line = String.trim (input_line ic) in
        if (not found_header) && line = header_start then
          loop true acc
        else if found_header && line <> "" then
          let a1, a2 =
            Scanf.sscanf line "%s\t\t%s\t\t%f"
              (fun allele1 allele2 conf ->
                { hla_class  = allele_to_hla_class allele1
                ; allele     = allele1
                ; qualifier = ""
                ; confidence = 1.0 -. conf
                },
                { hla_class  = allele_to_hla_class allele2
                ; allele     = allele2
                ; qualifier = ""
                ; confidence = 1.0 -. conf   (* Don't REALLY understand their metric *)
                })
          in
          let nacc =
            match equal_pairs with
            | `ReportAll -> a1 :: a2 :: acc
            | `FirstPair -> if acc = [] then a1 :: a2 :: [] else acc
            | `Unique    ->
              let add_if_new a set = if List.mem a ~set then set else a :: set in
              add_if_new a2 (add_if_new a1 acc)
          in
          loop true nacc
        else
          loop found_header acc
      with End_of_file ->
        close_in ic;
        acc
    in
    let alleles = loop false [] in
    run, alleles

  let scan_directory ?equal_pairs dir =
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
                  | None -> ms, ds
                  | Some m -> (df, m) :: ms, ds)
        in
        loop (matches @ acc) (dirs @ t)
    in
    loop [] [dir]
    |> List.map ~f:(parse ?equal_pairs)
    |> join_by_fst
    |> List.sort ~cmp:compare  (* sort by keys aka runs *)


end (* Athlates *)

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
