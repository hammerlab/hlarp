
open Std

module SMap = Map.Make (struct type t = string let compare = compare end)
module SSet = Set.Make (struct type t = string let compare = compare end)

let to_prefix name lst =
  if List.length lst > 1 then
    (fun i -> sprintf "%s%d" name (i + 1))
  else
    fun _ -> name

let nested_maps seqlst optlst athlst prolst filelst =
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
  let load_file file =
    let ic = open_in file in
    let _hdr = input_line ic in
    let rec loop m =
      try
        let line = input_line ic in
        let run, info =
          match Re.split (Re.char ',' |> Re.compile) line with
          | cls_s :: all_s :: qul_s :: con_s :: run_s :: [] ->
              run_s
              , { hla_class =
                    if cls_s = "1" then I else
                      if cls_s = "2" then II else
                        invalid_arg (sprintf "unrecognized cls: %s" cls_s)
                ; allele = all_s
                ; qualifier = qul_s
                ; confidence = float_of_string_nanable con_s
              }
          | _ -> invalid_arg
                    (sprintf "can't parse Hlarp input line %s from %s"
                      line file)
        in
        match SMap.find run m with
        | (f, lst) when f = file -> loop (SMap.add run (file, info :: lst) m)
        | exception Not_found    -> loop (SMap.add run (file, info :: []) m)
        | (nf, _)                ->
            invalid_arg (sprintf "encountered another file %s while parsing: %s" nf file)
      with End_of_file ->
        m
    in
    loop SMap.empty
  in
  let add_hlarp_files filelst run_map =
    List.fold_left filelst ~init:run_map ~f:(fun m file ->
      let fm = load_file file in
      let fmm = SMap.map (fun fl -> fl :: []) fm in
      SMap.union (fun _run l1 l2 -> Some (l1 @ l2)) m fmm)
  in
  SMap.empty
  |> add_to_run_map "seq2HLA"  Seq2HLA.scan_directory seqlst
  |> add_to_run_map "OptiType" OptiType.scan_directory optlst
  |> add_to_run_map "ATHLATES" Athlates.scan_directory athlst
  |> add_to_run_map "Prohlatype" Prohlatype.scan_directory prolst
  |> add_hlarp_files filelst
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

let jaccard ?(zero_on_empty=true) s1 s2 =
  let n1 = SSet.cardinal s1 in
  let n2 = SSet.cardinal s2 in
  if n1 = 0 && n2 = 0 then
    1.
  else if zero_on_empty && (n1 = 0 || n2 = 0) then
    0.
  else
    let inter = SSet.inter s1 s2 |> SSet.cardinal in
    let union = SSet.union s1 s2 |> SSet.cardinal in
    (float inter) /. (float union)

let to_filter = function
  | None                  -> fun _ -> true
  | Some (`HLAClass c)    -> fun ai -> ai.hla_class = c
  | Some (`LociPrefix p)  -> fun ai -> Re.execp (Re_posix.compile_pat ("^" ^ p)) ai.allele

let compute_mean_jaccard ?by ?(count_homozygous_2x=true) select typer_assoc =
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
  let sum_jaccard_index =
    fold_over_all_pairs as_sets ~init:0.
      ~f:(fun s ((_t1,s1), (_t2,s2)) -> s +. jaccard s1 s2)
  in
  let nf = float (List.length as_sets) in
  let num_pairs = nf *. (nf -. 1.) /. 2. in
  sum_jaccard_index /. num_pairs

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
    to ignoring the distinction.
    ?loci: Allow more than one HLA loci to do the analysis on, superseding
    any classes argument, but default to ignoring the distinction. Specify
    a string prefix that is used group alleles. *)
let output ?resolution ?classes ?loci oc nested_map_output =
  let select = select_allele ?resolution in
  let cmj, group, default_indent, suffix, nc =
    match loci with
    | Some loci_lst ->
        (fun typer_assoc ->
          List.map loci_lst ~f:(fun loci ->
            compute_mean_jaccard ~by:(`LociPrefix loci) select typer_assoc))
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
              (fun typer_assoc -> [ compute_mean_jaccard select typer_assoc ])
              , (fun typer_assoc -> [ None, group_similarities select typer_assoc ])
              , ""
              , "y"
              , 1
          | Some hla_classes ->
              (fun typer_assoc ->
                List.map hla_classes ~f:(fun hla_class ->
                  compute_mean_jaccard ~by:(`HLAClass hla_class) select typer_assoc))
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
            fprintf oc "%s\tjaccard similarit%s: %s\n" run suffix (float_lst_to_str jc_lst);
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
  fprintf oc "Average jaccard similarit%s across runs: %s\n" suffix (float_lst_to_str avgs)
