# Multiple Testing Across Genes and Contrasts

Function [`decideTests`](https://rdrr.io/pkg/limma/man/decideTests.html)
is not supported yet for a
[`PhyloMArrayLM`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md),
and will throw an error.

Function
[`classifyTestsF`](https://rdrr.io/pkg/limma/man/classifytestsF.html) is
not supported yet for a
[`PhyloMArrayLM`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md),
and will throw an error.

## Usage

``` r
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

  further parameters to be passed to
  [`classifyTestsF`](https://rdrr.io/pkg/limma/man/classifytestsF.html).
