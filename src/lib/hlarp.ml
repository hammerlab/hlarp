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
    let run = Re.Group.get grp 1 in
    let time_stamp = Re.Group.get grp 2 in
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
      (run, time_stamp), lst

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
        | (((first_run,_timestamp_dir),first_res) :: _) as lst ->
          match combine_multiday with
          | `TakeLast ->
              List.fold_left lst ~init:(first_run, first_res, [])
                ~f:(fun (previous_run, prev_res, acc) ((next_run, _timestamp_dir), lst) ->
                      if previous_run = next_run then
                        (next_run, lst, acc)  (* ignore the previous run's results. *)
                      else
                        (next_run, lst, (previous_run, prev_res) :: acc))
              |> (fun (p, l, acc) -> (p, l) :: acc))

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

module Compare = struct

  module SMap = Map.Make (struct type t = string let compare = compare end)
  module SSet = Set.Make (struct type t = string let compare = compare end)

  let to_prefix name lst =
    if List.length lst > 1 then
      (fun i -> sprintf "%s%d" name (i + 1))
    else
      fun _ -> name

  let nested_maps seqlst optlst athlst =
    let add_to_run_map name scan dirlst run_map =
      let prefix = to_prefix name dirlst in
      List.foldi dirlst ~init:run_map ~f:(fun i rmp dir ->
          let p = prefix i in
          scan dir
          |> List.fold_left ~init:rmp ~f:(fun rmp (run, info_lst) ->
              match SMap.find run rmp with
              | lst                 -> SMap.add run ((p,info_lst) :: lst) rmp
              | exception Not_found -> SMap.add run ((p,info_lst) :: []) rmp))
    in
    SMap.empty
    |> add_to_run_map "seq2HLA"  Seq2HLA.scan_directory seqlst
    |> add_to_run_map "OptiType" OptiType.scan_directory optlst
    |> add_to_run_map "ATHLATES" Athlates.scan_directory athlst
    |> SMap.bindings

  let fold_over_all_pairs ~f ~init lst =
    let rec loop init = function
      | []     -> init
      | h :: t -> loop (List.fold_left ~init ~f:(fun acc te -> f acc (h, te)) t) t
    in
    loop init lst

  let colon_regex = Re_posix.compile_pat ":"

  let select_allele ?resolution ai =
    match resolution with
    | None -> ai.allele
    | Some n -> List.take (Re.split colon_regex ai.allele) n
                |> String.concat ":"

  let jacard s1 s2 =
    if SSet.is_empty s1 && SSet.is_empty s2 then
      1.
    else
      let inter = SSet.inter s1 s2 |> SSet.cardinal in
      let union = SSet.union s1 s2 |> SSet.cardinal in
      (float inter) /. (float union)

  let to_filter = function
    | None                  -> fun _ -> true
    | Some (`HLAClass c)    -> fun ai -> ai.hla_class = c
    | Some (`LociPrefix p)  -> fun ai -> Re.execp (Re_posix.compile_pat ("^" ^ p)) ai.allele

  let compute_mean_jacard ?by ?(count_homozygous_2x=true) select typer_assoc =
    let filter = to_filter by in
    let set_of_ailst lst =
      List.fold_left lst ~init:SSet.empty
        ~f:(fun s ai ->
              if filter ai then
                let ai_sel = select ai in
                if count_homozygous_2x && SSet.mem ai_sel s then
                  SSet.add (ai_sel ^ "2") s
                else
                  SSet.add (select ai) s
              else
                s)
    in
    let as_sets =
      List.map ~f:(fun (typer, allele_lst) ->
        let s = set_of_ailst allele_lst in
        (typer, s)) typer_assoc
    in
    let sum_jacard_index =
      fold_over_all_pairs as_sets ~init:0.
        ~f:(fun s ((_t1,s1), (_t2,s2)) -> s +. jacard s1 s2)
    in
    let nf = float (List.length as_sets) in
    let num_pairs = nf *. (nf -. 1.) /. 2. in
    sum_jacard_index /. num_pairs

  let count_consecutive_doubles = function
    | []      -> []
    | h :: [] -> [ h, 1]
    | h :: t  ->
      let rec loop p c acc = function
        | []     -> List.rev ((p, c) :: acc)
        | h :: t when h = p -> loop p (c + 1) acc t
        | h :: t            -> loop h 1 ((p,c) :: acc) t
      in
      loop h 1 [] t

  let compress_counts =
    List.map ~f:(fun (a,c) -> if c < 2 then a else sprintf "%sx%d" a c)

  let group_similarities ?by select typer_assoc =
    let filter = to_filter by in
    List.fold_left typer_assoc ~init:SMap.empty
      ~f:(fun m (typer, allele_lst) ->
        List.fold_left allele_lst ~init:m ~f:(fun m ai ->
          if not (filter ai) then
            m
          else
            let ai_sel = select ai in
            match SMap.find ai_sel m with
            | lst                 -> SMap.add ai_sel (typer :: lst) m
            | exception Not_found -> SMap.add ai_sel (typer :: []) m))
    |> SMap.bindings
    |> List.sort ~cmp:compare
    |> List.map ~f:(fun (a, l) ->
        (a, count_consecutive_doubles l |> compress_counts))

  (* ?classes: Allow more than one HLA_class to do the analysis on, but default
      to ignoring the distinction. *)
  let output ?resolution ?classes ?loci oc nested_map_output =
    let select = select_allele ?resolution in
    let cmj, group, default_indent, suffix, nc =
      match loci with
      | Some loci_lst ->
          (fun typer_assoc ->
            List.map loci_lst ~f:(fun loci ->
              compute_mean_jacard ~by:(`LociPrefix loci) select typer_assoc))
          , (fun typer_assoc ->
            List.map loci_lst ~f:(fun loci ->
              (Some (loci ^ ":\t")
              , group_similarities ~by:(`LociPrefix loci) select typer_assoc)))
          , "\t"
          , "ies"
          , List.length loci_lst
      | None ->
          begin
            match classes with
            | None ->
                (fun typer_assoc -> [ compute_mean_jacard select typer_assoc ])
                , (fun typer_assoc -> [ None, group_similarities select typer_assoc ])
                , ""
                , "y"
                , 1
            | Some hla_classes ->
                (fun typer_assoc ->
                  List.map hla_classes ~f:(fun hla_class ->
                    compute_mean_jacard ~by:(`HLAClass hla_class) select typer_assoc))
                , (fun typer_assoc ->
                    List.map hla_classes ~f:(fun hla_class ->
                      (Some (hla_class_to_string hla_class ^ ":\t")
                      , group_similarities ~by:(`HLAClass hla_class) select typer_assoc)))
                , "\t"
                , "ies"
                , List.length hla_classes
          end
    in
    let float_lst_to_str l = String.concat " " (List.map ~f:(sprintf "%0.2f") l) in
    let ss, n =
      List.fold_left nested_map_output ~init:(List.init nc ~f:(fun _ -> 0.), 0)
        ~f:(fun (jc_s, jc_n) (run, typer_assoc) ->
              let jc_lst = cmj typer_assoc in
              fprintf oc "%s\tjacard similarit%s: %s\n" run suffix (float_lst_to_str jc_lst);
              let group = group typer_assoc in
              List.iter group ~f:(fun (cls_opt, by_cls_lst) ->
                List.iteri by_cls_lst ~f:(fun i (a, tlst) ->
                  let s1 =
                    if i = 0 then Option.value ~default:default_indent cls_opt
                    else default_indent
                  in
                  fprintf oc "\t%s%s\t%s\n" s1 a (String.concat ";" tlst)));
              List.map2 ~f:(+.) jc_s jc_lst, jc_n + 1)
    in
    let nf = float n in
    let avgs = List.map ~f:(fun x -> x /. nf) ss in
    fprintf oc "Average jacard similarit%s across runs: %s\n" suffix (float_lst_to_str avgs)

end (* Compare *)

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
