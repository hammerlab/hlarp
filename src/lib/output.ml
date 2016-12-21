
open Std

let nan_is_empty f = if f <> f then "" else sprintf "%f" f

let out_channel oc run_assoc =
  fprintf oc "class,allele,qualifier,confidence,typer specific,run\n";
  List.iter run_assoc ~f:(fun (run, lst) ->
    List.sort ~cmp:compare lst
    |> List.iter ~f:(fun {hla_class; allele; qualifier; confidence; typer_spec} ->
      fprintf oc "%s,%s,%s,%s,%s,%s\n"
        (hla_class_to_string hla_class)
        allele qualifier (nan_is_empty confidence) typer_spec
        run))

