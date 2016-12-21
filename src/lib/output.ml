
open Std

let out_channel oc run_assoc =
  fprintf oc "%s,run\n" Info.csv_header;
  List.iter run_assoc ~f:(fun (run, lst) ->
    List.sort ~cmp:compare lst
    |> List.iter ~f:(fun i ->
        Info.output_as_csv_row oc i;
        fprintf oc ",%s\n" run))
