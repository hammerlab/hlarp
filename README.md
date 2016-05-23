# hlarp

Normalize HLA typing output of [Seq2HLA](https://bitbucket.org/sebastian_boegel/seq2hla),
[Optitype](https://github.com/FRED-2/OptiType/) & [ATHLATES](https://www.broadinstitute.org/scientific-community/science/projects/viral-genomics/athlates).

`Hlarp` is also provided as a library.

Pronounced "heh-larp" (or with a soft 'h' : "larp").

## Installation

The easiest way to install the tool is to "[opam](http://opam.ocaml.org/) pin" this repository:

``` shell
opam pin add -k git https://github.com/hammerlab/hlarp
```

Now you can run hlarp with `hlarp`.

## Using hlarp

### Command line

#### Aggregating results

``` shell
hlarp seq2HLA /path/to/seq2HLA/results directory > results.csv
hlarp optitype /path/to/OptiType/results directory > results.csv
```

The output is a csv

| class | allele           | qualifier | confidence | run    |
|-------|------------------|-----------|------------|--------|
| 1     | A*02:05:01       |           | 0.500000   | 120013 |
| 1     | A*30:01:01       |           | 0.500000   | 120013 |
| 1     | B*57:03:01       |           | 1.000000   | 120013 |
| 1     | B*57:03:01       |           | 1.000000   | 120013 |
| 1     | C*18:02          |           | 1.000000   | 120013 |
| 1     | C*18:02          |           | 1.000000   | 120013 |
| 2     | DRB1*03:02:01    |           | 0.500000   | 120013 |
| 2     | DRB1*15:03:01:01 |           | 0.250000   | 120013 |
| 2     | DRB3*01:01:02:01 |           | 0.500000   | 120013 |
| 2     | DRB3*01:01:02:02 |           | 0.500000   | 120013 |
| 2     | DRB5*01:01:01    |           | 1.000000   | 120013 |
| 2     | DRB5*01:01:01    |           | 1.000000   | 120013 |
| 1     | A*31:01:02       |           | 0.500000   | 120021 |
| 1     | A*31:01:13       |           | 0.500000   | 120021 |
| 1     | C*02:10          |           | 1.000000   | 120021 |
| 1     | C*02:10          |           | 1.000000   | 120021 |
| 2     | DRB1*15:03:01:01 |           | 0.500000   | 120021 |
| 2     | DRB1*15:03:01:02 |           | 0.500000   | 120021 |
| 2     | DRB3*02:02:01:01 |           | 0.500000   | 120021 |
| 2     | DRB3*02:02:01:02 |           | 0.500000   | 120021 |
| 2     | DRB5*01:01:01    |           | 1.000000   | 120021 |
| 2     | DRB5*01:01:01    |           | 1.000000   | 120021 |
| 1     | A*68:02:01:01    |           | 0.166667   | 120074 |
| 1     | A*74:01          |           | 0.500000   | 120074 |
| 1     | B*15:03:01       |           | 0.500000   | 120074 |
| 1     | B*15:16:01       |           | 0.500000   | 120074 |
| 1     | C*02:10          |           | 0.500000   | 120074 |
| 1     | C*16:01:01       |           | 0.500000   | 120074 |
| 2     | DRB1*01:02:01    |           | 0.500000   | 120074 |
| 2     | DRB1*14:54:01    |           | 0.500000   | 120074 |


Some columns are left empty due to the nature of the HLA-typer.

#### Comparison

```shell
hlarp compare --resolution 2 -l A -l B -l C --loci DRB1 --optitype /path/to/optitype/results -a /path/to/ATHLATES/results/
```

Will generate this kind of report:

```
120013	jacard similarities: 1.00 1.00 0.00 0.00
	A:	A*02:05	OptiType;ATHLATES
		A*30:01	OptiType;ATHLATES
	B:	B*57:03	OptiTypex2;ATHLATESx2
	C:	C*18:01	OptiTypex2
		C*18:02	ATHLATESx2
	DRB1:	DRB1*03:02	ATHLATES
		DRB1*15:03	ATHLATES
120021	jacard similarities: 0.33 0.00 0.33 0.00
	A:	A*02:01	OptiType
		A*31:01	OptiType;ATHLATESx2
	B:	B*15:03	OptiType
		B*45:01	OptiType
	C:	C*02:10	OptiType;ATHLATESx2
		C*16:01	OptiType
	DRB1:	DRB1*15:03	ATHLATESx2
120074	jacard similarities: 1.00 1.00 1.00 0.00
	A:	A*68:02	OptiType;ATHLATES
		A*74:01	OptiType;ATHLATES
	B:	B*15:03	OptiType;ATHLATES
		B*15:16	OptiType;ATHLATES
	C:	C*02:10	OptiType;ATHLATES
		C*16:01	OptiType;ATHLATES
	DRB1:	DRB1*01:02	ATHLATES
		DRB1*14:54	ATHLATES
...
Average jacard similarities across runs: 0.54 0.60 0.51 0.07
```

### From OCaml top

You can also use the `hlarp` module:

``` ocaml
open Hlarp
```
