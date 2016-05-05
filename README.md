# hlarp

Normalize HLA typing output of Seq2HLA and Optitype.

`Hlarp` is also provided as a library.

Pronounced "heh-larp".

## Installation

``` shell
make setup # you don't need to do this more than once
make
```

Now you can run hlarp with `./hlarp_cli.native`.

## Using hlarp

``` shell
./hlarp_cli.native seq2HLA /path/to/seq2HLA results directory > results.csv
./hlarp_cli.native optitype /path/to/seq2HLA results directory > results.csv
```
