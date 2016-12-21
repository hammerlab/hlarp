
open Std

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
          Scanf.sscanf line "%s\t\t%s\t\t%s"
            (fun allele1 allele2 confidence ->
              { hla_class  = allele_to_hla_class allele1
              ; allele     = allele1
              ; qualifier  = ""
              ; confidence = float_of_string_nanable confidence
              ; typer_spec = ""
              },
              { hla_class  = allele_to_hla_class allele2
              ; allele     = allele2
              ; qualifier  = ""
              ; confidence = float_of_string_nanable confidence
              ; typer_spec = ""
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

