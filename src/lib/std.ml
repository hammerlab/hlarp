
include Nonstd
open Printf

module StringMap = Map.Make (struct type t = string let compare = compare end)

let join_by_fst lst =
  List.fold_left lst ~f:(fun acc (run, inf) ->
    match StringMap.find run acc with
    | exception Not_found -> StringMap.add run inf acc
    | clst -> StringMap.add run (clst @ inf) acc)
        ~init:StringMap.empty
  |> StringMap.bindings

let float_of_string_nanable s = try float_of_string s with Failure _ -> nan

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
