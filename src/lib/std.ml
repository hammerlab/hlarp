
include Nonstd
open Printf

let invalid_argf fmt = ksprintf invalid_arg fmt

module StringMap = Map.Make (struct type t = string let compare = compare end)

(* Combine assocation lists, with list values, by their string keys. *)
let join_by_fst lst =
  List.fold_left lst ~init:StringMap.empty ~f:(fun acc (r, inf) ->
      match StringMap.find r acc with
      | exception Not_found -> StringMap.add r inf acc
      | clst -> StringMap.add r (clst @ inf) acc)
  |> StringMap.bindings

let float_of_string_nanable s = try float_of_string s with Failure _ -> nan

type hla_class = | I | II

let hla_class_to_string = function | I -> "1" | II -> "2"

(* This type is a place holder in case we want to
   make a full ADT of MHC alleles in the future. *)
type allele = string

(** A sample (aka run) is an obvservation for which we may have multiple
    [info]'s. *)
type sample = string

module SampleMap = Map.Make (struct type t = sample let compare = compare end)
module AlleleMap = Map.Make (struct type t = allele let compare = compare end)


let nan_is_empty f = if f <> f then "" else sprintf "%f" f

module Info = struct

  type t =
    { hla_class   : hla_class
    ; allele      : allele
    ; qualifier   : string
    ; confidence  : float
    (* store typer specific information in free form text. *)
    ; typer_spec  : string
    }

  let to_string { hla_class; allele; qualifier; confidence; typer_spec} =
    sprintf "%s-%s-%s-%f-%s" (hla_class_to_string hla_class)
      allele qualifier confidence typer_spec

  let csv_header = "class,allele,qualifier,confidence,typer specific"

  let output_as_csv_row oc {hla_class; allele; qualifier; confidence; typer_spec} =
    fprintf oc "%s,%s,%s,%s,%s"
      (hla_class_to_string hla_class)
      allele
      qualifier
      (nan_is_empty confidence)
      typer_spec

  let comma_regex = Re.(compile (char ','))

  let input_as_csv_row ?(file="") ic =
    let line = input_line ic in
    match Re.split comma_regex line with
    | cls_s :: allele :: qualifier :: con_s :: typer_spec :: tl ->
        { hla_class =
          if cls_s = "1" then I else
            if cls_s = "2" then II else
              invalid_argf "unrecognized cls: %s" cls_s
        ; allele
        ; qualifier
        ; confidence  = float_of_string_nanable con_s
        ; typer_spec
        }, tl
    | _ ->
        invalid_argf "can't parse Hlarp input line %s from %s" line file

end (* Info *)

module InfoMap = Map.Make (struct type t = Info.t let compare = compare end)
