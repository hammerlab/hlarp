
open Std

let out_channel oc sample_assoc =
  fprintf oc "%s,sample\n" Info.csv_header;
  List.iter sample_assoc ~f:(fun (sample, lst) ->
    List.sort ~cmp:compare lst
    |> List.iter ~f:(fun i ->
        Info.output_as_csv_row oc i;
        fprintf oc ",%s\n" sample))
