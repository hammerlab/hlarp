
open Std

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
      { allele
      ; hla_class  = I
      ; qualifier  = ""
      ; confidence = nan
      ; typer_spec = ""
      }
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

