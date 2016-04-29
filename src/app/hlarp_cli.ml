(** A tool to normalize HLA typing output.

class,allele,qualifier,confidence
1,HLA-A*02:01,G,0.9
*)

open Nonstd
open Printf

open Hlarp

let project_name = "hlarp"

let () =
  let open Cmdliner in
  let version = "0.0.0" in
  let help_cmd =
    let doc = "HLA typing output normalizer." in
    let bug =
      sprintf "Browse and report new issues at <https://github.com/hammerlab/%s>."
        project_name
    in
    let man =
      [ `S "AUTHORS"
      ; `P "Leonid Rozenberg <leonidr@gmail.com>"
      ; `Noblank
      ; `S "BUGS"
      ; `P bug
      ]
    in
    Term.(ret (pure (`Help (`Plain, None)))
        , info project_name ~version ~doc ~man)
  in
  let directory_arg ~tool ~suffix =
    let doc =
      sprintf "Directory to find %s output.\
                      \n%s will search for files with %s suffix.\
                      Defaults to looking in the current dir."
              tool project_name suffix
    in
    Arg.(value & pos 0 dir "." & info ~doc [])
  in
  let seq2HLA =
    let tool = "seq2HLA" in
    Term.(const (fun dir ->
            Seq2HLA.scan_directory dir
            |> Output.out_channel stdout)
          $ (directory_arg ~tool ~suffix:Seq2HLA.suffix)
        , info tool ~doc:"Parse seq2HLA output")
  in
  let optitype =
    let tool = "optitype" in
    Term.(const (fun dir ->
            OptiType.scan_directory dir
            |> Output.out_channel stdout)
          $ (directory_arg ~tool ~suffix:OptiType.suffix)
        , info tool ~doc:"Parse OptiType output")
  in
  let cmds = [seq2HLA; optitype] in
  match Term.eval_choice help_cmd cmds with
  | `Ok () -> ()
  | `Error _ -> failwith "cmdliner error"
  | `Version | `Help -> exit 0
