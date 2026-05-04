# Methods for class PhyloMArrayLM

Methods for class PhyloMArrayLM

## Usage

``` r
getParameters(object, consensus = TRUE)

# S4 method for class 'PhyloMArrayLM'
getParameters(object, consensus = TRUE)

plotParameters(object, ...)

# S4 method for class 'PhyloMArrayLM'
plotParameters(object, ...)

# S4 method for class 'PhyloMArrayLM'
show(object)

logLikelihood(object)

# S4 method for class 'PhyloMArrayLM'
logLikelihood(object)

consensusTree(object)

# S4 method for class 'PhyloMArrayLM'
consensusTree(object)
```

## Arguments

- object:

  an object of class
  [`PhyloMArrayLM-class`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md)

- consensus:

  if `TRUE` (the default), `getParameters` returns the consensus
  parameters as a named vector. Otherwise, it returns the individual
  fits for all genes.

- ...:

  further parameters to be based to
  [`hist`](https://rdrr.io/r/graphics/hist.html), including `breaks`.

## Value

Depending on the method, a named vector, a data frame, or a plot.

## See also

[`PhyloMArrayLM-class`](https://pbastide.github.io/phyloDE/reference/PhyloMArrayLM-class.md),
[`phylolmFit`](https://pbastide.github.io/phyloDE/reference/phylolmFit.md)
