# Compute 1 - rho

rhoprime = 1 - rho = (1 - exp(-2\*alpha\*t_H)) / (2 \* alpha \* t_H) is
the variance of the OU over the variance of the BM. Taken from
Cornuault, 2023, Syst. Biol. It is the fraction of the variance that can
be explained by "neutral" BM evolution. When alpha goes to 0, rhoprime
goes to 1: trait can be explained by "neutral" BM process. When alpha
goes to Inf, rhoprime goes to 0: trait is deterministic.

## Usage

``` r
rho_prime(alpha, t_tree, tol = .Machine$double.eps)
```

## Arguments

- alpha:

  the selection strength of the process

- t_tree:

  the total height of the tree

- tol:

  tolerence value for calling alpha = 0.

## Value

Value of rhoprime
