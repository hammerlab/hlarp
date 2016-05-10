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

let hla_class_to_string = function | I -> "1" | II -> "2"

(* This type is a place holder in case we want to
   make a full ADT of MHC alleles in the future. *)
type allele = string

type info =
  { hla_class   : hla_class
  ; allele      : allele
  ; qualifier   : string
  ; confidence  : float
  }

let info_to_string { hla_class; allele; qualifier; confidence} =
  sprintf "%s-%s-%s-%f" (hla_class_to_string hla_class)
    allele qualifier confidence

module InfoMap = Map.Make (struct type t = info let compare = compare end)

module Seq2HLA = struct

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

  let parse ?(equal_pairs=`MostLikelyPair) (fname, re_group) =
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
              (fun allele1 allele2 confidence ->
                { hla_class  = allele_to_hla_class allele1
                ; allele     = allele1
                ; qualifier  = ""
                ; confidence
                },
                { hla_class  = allele_to_hla_class allele2
                ; allele     = allele2
                ; qualifier  = ""
                ; confidence
                })
          in
          loop true (a1 :: a2 :: acc)
        else
          loop found_header acc
      with End_of_file ->
        close_in ic;
        acc
    in
    let alleles = loop false [] in
    let alleles =
      match equal_pairs with
      | `ReportAll      -> alleles
                        (* Will fail silently if < 2, I think ok. *)
      | `FirstPair      -> List.take (List.rev alleles) 2
      | `Unique         -> List.dedup alleles
      | `MostLikelyPair ->
          let totalf = float (List.length alleles) in
          List.fold_left alleles ~f:(fun m a ->
            match InfoMap.find a m with
            | exception Not_found -> InfoMap.add a 1 m
            | n -> InfoMap.add a (n + 1) m) ~init:InfoMap.empty
          |> InfoMap.bindings
          |> List.sort ~cmp:(fun (_a1, c1) (_a2, c2) -> compare c2 c1 (* reverse *))
          |> (fun lst -> List.take lst 2)
          |> List.map ~f:(fun (a, c) -> { a with confidence = (float c) /. totalf })
          |> function
              | a :: [] -> [a; a] (* Preserve pair to signal homozygosity *)
              | lst     -> lst
    in
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

module CompareByLocus = struct

  type locus =
    | A
    | B
    | C
    | DRB1
    | DRB3
    | DRB4
    | DRB5

  let locus_to_string = function
    | A     -> "A"
    | B     -> "B"
    | C     -> "C"
    | DRB1  -> "DRB1"
    | DRB3  -> "DRB3"
    | DRB4  -> "DRB4"
    | DRB5  -> "DRB5"

  let known_loci = [ A; B; C; DRB1; DRB3; DRB4; DRB5 ]

  let allele_to_locus allele =
    try
      let i = String.index allele '*' in
      match String.sub allele 0 i with
      | "A"    -> Some A
      | "B"    -> Some B
      | "C"    -> Some C
      | "DRB1" -> Some DRB1
      | "DRB3" -> Some DRB3
      | "DRB4" -> Some DRB4
      | "DRB5" -> Some DRB5
      | _ -> eprintf "Unrecognized allele to locus conversion: %s" allele;
             None
    with Not_found ->
      eprintf "Unrecognized allele to locus conversion: %s" allele;
      None

  (* An ad-hoc grouping operation. *)
  let group by select lst =
    let rec loop acc = function
      | []     -> acc
      | h :: t ->
        let b = by h in
        let s = select h in
        let like, not_like =
          List.partition_map t ~f:(fun l ->
            if by l = b then `Fst (select l) else `Snd l)
        in
        (* The sort here allows us to compare lists *)
        let s_like = List.sort ~cmp:compare (s :: like) in
        let n = List.length s_like in
        loop ((n, b, s_like) :: acc) not_like
  in
  loop [] lst

  let by_loci lst =
    let by (_key, ai) = ai.allele in
    let select (key, _ai) = key in
    let lm =
      List.fold_left lst ~init:[] ~f:(fun lm (ai, k) ->
        match allele_to_locus ai.allele with
        | None   -> lm  (* ignore, warning triggered by failed conversion. *)
        | Some l -> (l, (k, ai)) :: lm)
    in
    let g1, _g2 =
      List.fold_left known_loci ~init:([], lm) ~f:(fun (acc, lm) locus ->
          let from_locus, not_from_locus =
            List.partition ~f:(fun (l,_) -> l = locus) lm
          in
          if from_locus = [] then
            (acc, not_from_locus) (* locus not report *)
          else
            let f1 = List.map ~f:snd from_locus in
            let f2 = group by select f1 in
            let in_locus = List.sort ~cmp:(fun (n1, _, _) (n2, _, _) -> compare n2 n1) f2 in
            (locus, in_locus) :: acc, not_from_locus)
    in
    List.rev g1
    (*|> snd
    |> List.rev  get results back in 'known_loci' order *)

  let output oc lst =
    List.iter lst ~f:(fun (locus, in_locus) ->
      Printf.fprintf oc "%s:\n" (locus_to_string locus);
      List.iter in_locus ~f:(fun (_n, allele, keys_lst) ->
        Printf.fprintf oc "\t%s:\t%s\n" allele (String.concat "," keys_lst)))
    (*  List.iter in_locus ~f:(fun key ->
        Printf.fprintf oc "\t%s:\t%s\n" allele (String.concat "," keys))) *)

end (* CompareByLocus *)

module Output = struct

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
