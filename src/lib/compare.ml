
open Std

module SampleMap = Map.Make (struct type t = sample let compare = compare end)
module AlleleMap = Map.Make (struct type t = allele let compare = compare end)

type directory = string
type scanner = directory -> (sample * Info.t list) list

type source = string

(* When we want to compare multiple runs from the same source type
   (ex. Seq2HLA dir1, dir2 ..) use this method to create a prefix function
   that appends a number to the source if there are multiple and
   create source. *)
let to_source (name : string) (lst : 'a list) : (int -> source) =
  if List.length lst > 1 then
    (fun i -> sprintf "%s%d" name (i + 1))
  else
    fun _ -> name

let add_to_sample_map ~source_name (scan : scanner) (dirlst : directory list) sample_map =
  let prefix = to_source source_name dirlst in
  List.foldi dirlst ~init:sample_map ~f:(fun i rmp dir ->
      let p = prefix i in
      scan dir
      |> List.fold_left ~init:rmp ~f:(fun rmp (sample, info_lst) ->
          match SampleMap.find sample rmp with
          | lst                 -> SampleMap.add sample ((p,info_lst) :: lst) rmp
          | exception Not_found -> SampleMap.add sample ((p,info_lst) :: []) rmp))

(* Use previously written Hlarp output files as comparison source.  *)
module HlarpOutputFiles = struct

  let load_file ~use_basename (basename, file) =
    let ic = open_in file in
    let _hdr = input_line ic in
    let rec loop m =
      try
        let info, tl = Info.input_as_csv_row ~file ic in
        let sample =
          match tl with
          | [s] -> s
          | []  -> invalid_argf "Missing sample name in: %s" file
          | lst -> invalid_argf "Too many sample names after info: %s, in %s."
                    (String.concat "," lst) file
        in
        let key = if use_basename then basename else file in
        match SampleMap.find sample m with
        | (f, lst) when f = key -> loop (SampleMap.add sample (key, info :: lst) m)
        | exception Not_found   -> loop (SampleMap.add sample (key, info :: []) m)
        | (nf, _)               -> invalid_arg
              (sprintf "encountered another file %s while parsing: %s" nf key)
      with End_of_file ->
        m
    in
    loop SampleMap.empty

  let setup filelst =
    let with_basenames =
      List.map filelst ~f:(fun file -> Filename.basename file, file)
      |> List.sort ~cmp:compare
    in
    let compare_by_first (b1, _) (b2, _) = compare b1 b2 in
    let use_basename = not (List.contains_dup with_basenames ~compare:compare_by_first) in
    use_basename, with_basenames


  let add filelst sample_map =
    let use_basename, with_basenames = setup filelst in
    List.fold_left with_basenames ~init:sample_map ~f:(fun m bfp ->
      let fm = load_file ~use_basename bfp in
      let fmm = SampleMap.map (fun fl -> fl :: []) fm in
      SampleMap.union (fun _sample l1 l2 -> Some (l1 @ l2)) m fmm)

end

let hlarp_file_arg = "hlarp-file"
let seq2hl_file_arg = "seq2HLA"
let optitype_file_arg = "optitype"
let athlates_file_arg = "athlates"
let prohlatype_file_arg = "prohlatype"

let seq2HLA_source = "seq2HLA"
let optitype_source = "OptiType"
let athlates_source = "ATHLATES"
let prohlatype_source = "Prohlatype"

let argument_order ~use_basename =
  let rec loop s o a p acc = function
    | []
    | _ :: []       ->  List.rev acc
    | h :: t :: tl  ->
        if h = "--" ^ hlarp_file_arg then
          loop s o a p ((if use_basename then Filename.basename t else t) :: acc) tl
        else if h = "--" ^ seq2hl_file_arg then
          let source_name = if s = 0 then seq2HLA_source else seq2HLA_source ^ (string_of_int s) in
          loop (s + 1) o a p (source_name :: acc) tl
        else if h = "--" ^ optitype_file_arg then
          let source_name = if s = 0 then optitype_source else optitype_source ^ (string_of_int o) in
          loop s (o + 1) a p (source_name :: acc) tl
        else if h = "--" ^ athlates_file_arg then
          let source_name = if s = 0 then athlates_source else athlates_source ^ (string_of_int o) in
          loop s o (a + 1) p (source_name :: acc) tl
        else if h = "--" ^ prohlatype_file_arg then
          let source_name = if s = 0 then prohlatype_source else prohlatype_source ^ (string_of_int o) in
          loop s o a (p + 1) (source_name :: acc) tl
        else
          loop s o a p acc (t :: tl)
  in
  loop 0 0 0 0 [] (Array.to_list Sys.argv)

let nested_maps ?(seqlst=[]) ?(optlst=[]) ?(athlst=[]) ?(prolst=[]) ?(filelst=[]) () =
  let loaded_list =
    SampleMap.empty
    |> add_to_sample_map ~source_name:seq2HLA_source  Seq2HLA.scan_directory seqlst
    |> add_to_sample_map ~source_name:optitype_source OptiType.scan_directory optlst
    |> add_to_sample_map ~source_name:athlates_source Athlates.scan_directory athlst
    |> add_to_sample_map ~source_name:prohlatype_source Prohlatype.scan_directory prolst
    |> HlarpOutputFiles.add filelst
    |> SampleMap.bindings
  in
  let use_basename, _ = HlarpOutputFiles.setup filelst in
  let arg_order = argument_order ~use_basename in
  Printf.printf "Args: %s\n" (Sys.argv |> Array.to_list |> String.concat "; ");
  Printf.printf "Arg order: %s\n" (String.concat "; " arg_order);

  let find_in_arg source =
    let ns = String.length source in
    List.findi arg_order ~f:(fun _i s ->
        String.sub s 0 (min (String.length s) ns) = source)
  in
  List.map loaded_list ~f:(fun (k, assoc) ->
    k, List.sort assoc ~cmp:(fun (s1, _) (s2, _) ->
        compare (find_in_arg s1) (find_in_arg s2)))

let list_fold_over_all_pairs ~f ~init lst =
  let rec loop init = function
    | []     -> init
    | h :: t -> loop (List.fold_left ~init ~f:(fun acc te -> f acc h te) t) t
  in
  loop init lst

let colon_regex = Re_posix.compile_pat ":"

let to_allele ?resolution ai =
  match resolution with
  | None   -> ai.Info.allele
  | Some n -> List.take (Re.split colon_regex ai.Info.allele) n
              |> String.concat ":"

module AlleleSet = Set.Make (struct type t = allele let compare = compare end)

let jaccard ?(zero_on_empty=true) s1 s2 =
  let n1 = AlleleSet.cardinal s1 in
  let n2 = AlleleSet.cardinal s2 in
  if n1 = 0 && n2 = 0 then
    1.
  else if zero_on_empty && (n1 = 0 || n2 = 0) then
    0.
  else
    let inter = AlleleSet.inter s1 s2 |> AlleleSet.cardinal in
    let union = AlleleSet.union s1 s2 |> AlleleSet.cardinal in
    (float inter) /. (float union)

let loci_prefix_filter p ai =
  Re.execp (Re_posix.compile_pat ("^" ^ p)) ai.Info.allele

let hla_class_filter c ai =
  ai.Info.hla_class = c

let count_consecutive_doubles = function
  | []      -> []
  | h :: [] -> [ h, 1]
  | h :: t  ->
    let rec loop p c acc = function
      | []                -> List.rev ((p, c) :: acc)
      | h :: t when h = p -> loop p (c + 1) acc t
      | h :: t            -> loop h 1 ((p,c) :: acc) t
    in
    loop h 1 [] t

let compress_counts =
  List.map ~f:(fun (a,c) -> if c < 2 then a else sprintf "%sx%d" a c)

let loci_class_filters loci classes =
  match loci, classes with
  | Some loci, _              -> List.map loci ~f:(loci_prefix_filter)
  | None,      None           -> [ fun _ -> true ]
  | None,      (Some classes) -> List.map classes ~f:(hla_class_filter)

let specific_grouped_view ?loci ?classes ~typer1 ~ilist1 ~typer2 ~ilist2 =
  let filters = loci_class_filters loci classes in
  let pass_filters aiv = List.exists filters ~f:(fun f -> f aiv) in
  let to_fold typer =
    fun m (aik, aiv) ->
      if not (pass_filters aiv) then m else
        match AlleleMap.find aik m with
        | lst                 -> AlleleMap.add aik (typer :: lst) m
        | exception Not_found -> AlleleMap.add aik (typer :: []) m
  in
  List.fold_left ilist1 ~init:AlleleMap.empty ~f:(to_fold typer1)
  |> fun m -> List.fold_left ilist2 ~init:m ~f:(to_fold typer2)
  |> AlleleMap.bindings
  |> List.map ~f:(fun (a, l) ->
      (a, count_consecutive_doubles l |> compress_counts))
  |> List.sort ~cmp:(fun (_a1, s1) (_a2, s2) ->
      -1 * compare (List.length s1) (List.length s2))

let against_all_pairs eval nested_map_bindings =
  List.map nested_map_bindings ~f:(fun (sample, type_assoc) ->
    sample
    , list_fold_over_all_pairs type_assoc ~init:[]
        ~f:(fun acc ((typer1,_) as p1) ((typer2, _) as p2) ->
              let res = sprintf "%s vs %s" typer1 typer2 in
              let e   = eval p1 p2 in
              (res, e) :: acc))

type comparison =
  { metrics_eval  : float list
  ; grouped_view  : (allele * string list) list
  }

let args_to_projections ?loci ?classes mlst spec_gv data =
  let projects = loci_class_filters loci classes in
  let eval ((_t1, ilist1) as itassoc1) ((_t2, ilist2) as itassoc2) =
    { metrics_eval=
        List.map projects ~f:(fun fltr ->
          let f1 = List.filter ilist1 ~f:fltr in
          let f2 = List.filter ilist2 ~f:fltr in
          List.map mlst ~f:(fun k -> k f1 f2))
        |> List.concat
    ; grouped_view = spec_gv ~itassoc1 ~itassoc2
    }
  in
  against_all_pairs eval data

let remove_and_assoc el list =
  let rec loop acc = function
    | []                             -> raise Not_found
    | ((e, _) as p) :: t when e = el -> p, (List.rev acc @ t)
    | h :: t                         -> loop (h :: acc) t
  in
  loop [] list

(* Reduce the [info]'s to have a unique allele resolution.
   ex. A*01:01:01:01 -> A*01:01 and 2nd pair will get "x2" appended. *)
let keyed_by_allele_info_assoc ?resolution ?(label_homozygous=true) =
  let open Info in
  let to_allele = to_allele ?resolution in
  List.fold_left ~init:[] ~f:(fun acc ai ->
    let k = to_allele ai in
    if List.mem_assoc k acc then begin
      if label_homozygous then
        (k ^ "x2", ai) :: acc
      else
        let (_, aio), without = remove_and_assoc k acc in
        (k, {aio with confidence = aio.confidence +. ai.confidence }) :: without
    end else
      (k, ai) :: acc)

let set_of_assoc_keys l =
  (AlleleSet.of_list (List.map ~f:fst l))

let sum_confidence p =
  let open Oml in
  List.map p ~f:(fun (_, ai) -> ai.Info.confidence)
  |> Array.of_list
  |> Util.Array.sumf

let sum_snd p =
  let open Oml in
  List.map p ~f:snd
  |> Array.of_list
  |> Util.Array.sumf

let normalize sum assoc =
  let open Oml in
  if Util.is_degenerate sum || sum = 0. then
    let oon = 1. /. (float (List.length assoc)) in
    List.map ~f:(fun (k, _) -> k, oon) assoc
  else
    List.map ~f:(fun (k, ai) -> k, ai.Info.confidence /. sum) assoc

let describe_distr p =
  String.concat ";" (List.map p ~f:(fun (s, c) -> sprintf "\"%s\",%0.20f" s c))

let to_metric = function
  | `Jaccard ->
      begin fun asc1 asc2 ->
        jaccard (set_of_assoc_keys asc1) (set_of_assoc_keys asc2)
      end
      , true
  | `KLDiv   ->
      let open Oml in
      begin fun asc1 asc2 ->
        let p = normalize (sum_confidence asc1) asc1 in
        let q = normalize (sum_confidence asc2) asc2 in
        try
          let i = Statistics.Measures.discrete_kl_divergence ~d:1e-10 ~p ~q () in
          if i = infinity then begin
            let missing =
              List.filter_map p ~f:(fun (k, _) ->
                  if not (List.mem_assoc k q) then Some k else None)
            in
            eprintf "infinity! alleles found in p but not in q: %s \n\tp: %0.16f %s\n\tq: %0.16f %s\n"
              (String.concat "; " missing)
              (sum_snd p) (describe_distr p)
              (sum_snd q) (describe_distr q)
          end;
          i
        with (Invalid_argument m) ->
          eprintf "%s returning infinity. p: %0.20f %s vs q: %0.20f %s\n" m
            (sum_snd p) (describe_distr p)
            (sum_snd q) (describe_distr q);
          infinity
      end
      , false

(* ?classes: Allow more than one HLA_class to do the analysis on, but default
    to ignoring the distinction.
    ?loci: Allow more than one HLA loci to do the analysis on, superseding
    any classes argument, but default to ignoring the distinction. Specify
    a string prefix that is used group alleles.  *)
let analyze_samples ?loci ?classes ?resolution ?label_homozygous ?(metrics=[`Jaccard])
  nested_map_bindings =
    let ilist_map = keyed_by_allele_info_assoc ?resolution in
    let metricsf =
      List.map metrics ~f:(fun m ->
          let mtr, label_homozygous = to_metric m in
          fun alst1 alst2 ->
            let ilist1 = ilist_map ~label_homozygous alst1 in
            let ilist2 = ilist_map ~label_homozygous alst2 in
            mtr ilist1 ilist2)
    in
    let spgv ~itassoc1:(typer1, ilist1) ~itassoc2:(typer2, ilist2) =
      specific_grouped_view ?loci ?classes
        ~typer1 ~ilist1:(ilist_map ?label_homozygous ilist1)
        ~typer2 ~ilist2:(ilist_map ?label_homozygous ilist2)
    in
    args_to_projections ?loci ?classes metricsf spgv nested_map_bindings

let float_lst_to_str ~sep l =
  String.concat sep (List.map ~f:(sprintf "%0.2f") l)

let output_aggregates ?(summary_by=`Mean) oc smpl_width widths analyze_samples_output =
  let sep = "\t" in
  let _fst, fst_assoc = List.hd_exn analyze_samples_output in
  let num_comparisons = List.length fst_assoc in
  let num_metrics     = List.length (List.hd_exn fst_assoc |> snd).metrics_eval in
  let init = List.init num_comparisons ~f:(fun _ -> List.init num_metrics ~f:(fun _ -> [])) in
  let all_metrics =
    List.fold_left analyze_samples_output ~init ~f:(fun lst (_smpl, asc) ->
      List.map2 lst asc ~f:(fun acc (_, {metrics_eval;_}) ->
        List.map2 acc metrics_eval ~f:(fun lst v -> v :: lst)))
  in
  let open Oml.Statistics.Descriptive in
  fprintf oc "%-*s" smpl_width
    (match summary_by with
      | `Mean -> "mean"
      | `Median -> "median");
  List.map2 widths all_metrics ~f:(fun w lst ->
    let reduced =
      List.map lst ~f:(fun mtrlst ->
        let arr = Array.of_list mtrlst in
        match summary_by with
        | `Mean -> mean arr
        | `Median -> median arr) in
    sprintf "%-*s" w (float_lst_to_str ~sep:", " reduced))
  |> String.concat sep
  |> fprintf oc "%s\n"

let output ?(max_allele_rows_to_print=max_int) ?summary_by oc analyze_samples_output =
  let fst_smpl, fst_assoc = List.hd_exn analyze_samples_output in
  if fst_assoc = [] then fprintf oc "No comparisons\n" else
    let join_allele_str (a, lst) =
      sprintf "%s: %s" a (String.concat "," lst)
    in
    let widths =
      List.map ~f:(fun (s, r) ->
        List.fold_left r.grouped_view ~init:(0, String.length s)
          ~f:(fun (c, x) al ->
                if c > max_allele_rows_to_print then
                  (c + 1, x)
                else
                  (c + 1, max x (String.length (join_allele_str al))))
        |> snd) fst_assoc
    in
    let smpl_width = max (String.length "sample") (String.length fst_smpl + 2) in
    let sep = "\t" in
    let row_printer comp_assoc =
      List.map2 ~f:(fun w (s, _) -> sprintf "%*s" w s) widths comp_assoc
      |> String.concat sep
      |> fprintf oc "%s\n"
    in
    let value_row_printer comp_assoc =
      List.map2 ~f:(fun w (_, m) ->
        sprintf "%-*s" w (float_lst_to_str ~sep:", " m.metrics_eval)) widths comp_assoc
      |> String.concat sep
      |> fprintf oc "%s\n"
    in
    let setup_sub_row =
      List.map ~f:(fun (_, m) -> List.length m.grouped_view, m.grouped_view)
    in
    let sub_row_printer ssr row =
    List.map2 ~f:(fun w (ng, gv) ->
        if row < ng then
          sprintf "%*s" w (join_allele_str (List.nth_exn gv row))
        else
          sprintf "%*s" w " ") widths ssr
      |> String.concat sep
      |> fprintf oc "%s\n"
    in
    fprintf oc "%-*s" smpl_width "sample" ;
    row_printer fst_assoc;
    List.iter analyze_samples_output ~f:(fun (sample, comp_assoc) ->
      fprintf oc "%-*s " smpl_width sample;
      value_row_printer comp_assoc;
      let sub_rows = setup_sub_row comp_assoc in
      let mx_rows = List.fold_left ~init:min_int ~f:(fun x (r, _) -> max x r) sub_rows in
      let num_rows = min max_allele_rows_to_print mx_rows - 1 in
      for r = 0 to num_rows do
        fprintf oc "%*s " smpl_width " ";
        sub_row_printer sub_rows r
      done);
    output_aggregates oc smpl_width widths ?summary_by analyze_samples_output
