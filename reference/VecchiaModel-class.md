# Univariate Vecchia PrestoGP model class

This class is used to create univariate models with a likelihood
function conditioned on a subset of the observations (i.e., Vecchia
models). See
[`PrestoGPModel`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md)
for a description of the slots.

## See also

[`PrestoGPModel`](https://niehs.github.io/PrestoGP/reference/PrestoGPModel-class.md)

## Examples

``` r
pgp.model <- new("VecchiaModel", n_neighbors = 25)
```
