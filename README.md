# hlarp

Normalize HLA typing output of Seq2HLA and Optitype.

`Hlarp` is also provided as a library.

Pronounced "heh-larp".

## Installation

You can get the tool from oasis with `opam install hlarp`, or download it and
install from master locally, from the `hlarp` git directory, with the following:

``` shell
opam pin add -k git hlarp .
```

Now you can run hlarp with `hlarp`.

## Using hlarp

``` shell
hlarp seq2HLA /path/to/seq2HLA results directory > results.csv
hlarp optitype /path/to/seq2HLA results directory > results.csv
```

You can also use the `hlarp` module: 

``` ocaml
open Hlarp
...
```

