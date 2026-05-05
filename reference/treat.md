# Not Implemented Functions

Functions [`treat`](https://rdrr.io/pkg/limma/man/ebayes.html),
[`decideTests`](https://rdrr.io/pkg/limma/man/decideTests.html),
[`classifyTestsF`](https://rdrr.io/pkg/limma/man/classifytestsF.html)
are not supported yet for a
[`PhyloMArrayLM`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md),
and will throw an error.

## Usage

``` r
treat(fit, ...)

decideTests(fit, ...)

classifyTestsF(fit, ...)
```

## Arguments

- fit:

  a
  [`PhyloMArrayLM`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md)
  object, fitted using
  [`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md).

- ...:

  further parameters.
