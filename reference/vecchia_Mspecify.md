# Specify a multivariate Vecchia approximation

Specifies a multivariate Vecchia approximation for later use in
likelihood evaluation or prediction. This function does not depend on
parameter values, and only has to be run once before repeated likelihood
evaluations. This function is a multivariate version of
[`vecchia_specify`](https://rdrr.io/pkg/GPvecchia/man/vecchia_specify.html).

## Usage

``` r
vecchia_Mspecify(
  locs.list,
  m,
  locs.list.pred = NULL,
  dist.func = NULL,
  ordering.pred = c("obspred", "general"),
  pred.cond = c("independent", "general")
)
```

## Arguments

- locs.list:

  List of observed locations. Each each element should be a matrix
  containing the locs for the corresponding outcome variable.

- m:

  Number of nearby points to condition on.

- locs.list.pred:

  List of locations at which to make predictions. Each element should be
  a matrix containing the locs for the corresponding outcome variable.

- dist.func:

  Any distance function with a signature of dist(query_location,
  locations_matrix). Defaults to Euclidean distance.

- ordering.pred:

  Should "obspred" or "general" ordering be used for prediction? See
  [`vecchia_specify`](https://rdrr.io/pkg/GPvecchia/man/vecchia_specify.html).
  Defaults to "obspred".

- pred.cond:

  Should prediction conditioning be "general" or "independent"? See
  [`vecchia_specify`](https://rdrr.io/pkg/GPvecchia/man/vecchia_specify.html).
  Defaults to "independent".

## Value

An object that specifies the multivariate Vecchia approximation for
later use in likelihood evaluation or prediction.

## Details

This function should produce identical results to
[`vecchia_specify`](https://rdrr.io/pkg/GPvecchia/man/vecchia_specify.html)
for univariate problems, although it has fewer options. We recommend
that
[`vecchia_specify`](https://rdrr.io/pkg/GPvecchia/man/vecchia_specify.html)
be used in the univariate case.

## References

- Katzfuss, M., and Guinness, J. "A general framework for Vecchia
  approximations of Gaussian processes", Statistical Science (2021)
  36(1):124-141.

## See also

[`vecchia_specify`](https://rdrr.io/pkg/GPvecchia/man/vecchia_specify.html)

## Examples

``` r
data(soil)
soil <- soil[!is.na(soil[,5]),] # remove rows with NA's
locs <- as.matrix(soil[,1:2])
locsm <- list()
locsm[[1]] <- locsm[[2]] <- locs
soil.va <- vecchia_Mspecify(locsm, m=10)
```
