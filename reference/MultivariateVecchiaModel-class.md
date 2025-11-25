# Multivariate Vecchia PrestoGP model class

This class is used to create multivariate models with a likelihood
function conditioned on a subset of the observations (i.e., Vecchia
models). See
[`PrestoGPModel`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md)
for a description of the slots.

## See also

[`PrestoGPModel`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md)

## Examples

``` r
pgp.mmodel <- new("MultivariateVecchiaModel", n_neighbors = 25)
```
