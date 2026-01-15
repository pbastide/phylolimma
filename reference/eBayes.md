# Empirical Bayes Statistics for Differential Expression

Apply [`eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html) to the
result of function
[`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md).

Function [`treat`](https://rdrr.io/pkg/limma/man/ebayes.html) is not
supported yet for a
[`PhyloMArrayLM`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md),
and will throw an error.

## Usage

``` r
eBayes(fit, ...)

treat(fit, ...)
```

## Arguments

- fit:

  a
  [`PhyloMArrayLM`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md)
  object, fitted using
  [`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md).

- ...:

  further parameters to be passed to
  [`eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html) or
  [`treat`](https://rdrr.io/pkg/limma/man/ebayes.html).
