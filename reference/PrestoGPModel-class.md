# PrestoGPModel superclass

This is the superclass for PrestoGP models. All other types of PrestoGP
model classes (e.g.
[`VecchiaModel-class`](https://niehs.github.io/PrestoGP/reference/VecchiaModel-class.md)
and
[`FullModel-class`](https://niehs.github.io/PrestoGP/reference/FullModel-class.md)
are inherited from this class. Normally users should not create objects
of this class. Instead, they should use the appropriate inherited class
for the type of model they are fitting.

## Slots

- `covparams`:

  A numeric vector containing the parameters for the Matern model.

- `logparams`:

  A numeric vector containing the transformed versions of the Matern
  parameters (used internally for likelihood calculations).

- `beta`:

  A column matrix containing the regression coefficients.

- `vecchia_approx`:

  The output of the Vecchia specify function. See
  [`vecchia_specify`](https://rdrr.io/pkg/GPvecchia/man/vecchia_specify.html)
  and
  [`vecchia_Mspecify`](https://niehs.github.io/PrestoGP/reference/vecchia_Mspecify.md).

- `X_train`:

  A matrix containing the original predictors. This will be a "super
  matrix" for multivariate models. See
  [`superMatrix`](https://rdrr.io/pkg/psych/man/super.matrix.html).

- `Y_train`:

  A column matrix containing the original response values. After fitting
  the model, missing values will be replaced with the imputed values if
  impute.y is TRUE. See `link{prestogp_fit}`.

- `X_ndx`:

  A vector used to find the elements of beta corresponding to specific
  outcomes. The *i*th element of X_ndx is the index of the last element
  of beta corresponding to predictors for outcome *i*.

- `Y_bar`:

  A vector containing the means of each outcome.

- `Y_obs`:

  A logical vector used to track which values of Y are non-missing.

- `X_tilde`:

  The matrix of transformed predictors.

- `y_tilde`:

  The column matrix containing the transformed response values.

- `res`:

  A numeric vector of the residuals.

- `locs_train`:

  A list containing the location coordinates. Each element of the list
  corresponds to a different outcome. (The list will have length 1 for
  univariate models.)

- `penalty`:

  The type of penalized regression used. Should be one of "lasso",
  "relaxed", "MCP", or "SCAD".

- `linear_model`:

  The glmnet or ncvreg model. See
  [`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html),
  [`ncvreg`](https://pbreheny.github.io/ncvreg/reference/ncvreg.html),
  [`cv.glmnet`](https://glmnet.stanford.edu/reference/cv.glmnet.html),
  and
  [`cv.ncvreg`](https://pbreheny.github.io/ncvreg/reference/cv.ncvreg.html).

- `converged`:

  Did the model fitting process converge (boolean)?

- `LL_Vecchia_krig`:

  The value of the negative log likelihood function after optimization.

- `error`:

  Penalized model error. See References for details.

- `n_neighbors`:

  Number of neighbors to condition on for the Vecchia approximation.
  Ignored for full models.

- `min_m`:

  Minimum permissible number of neighbors.

- `alpha`:

  Parameter alpha for glmnet or ncvreg. See
  [`glmnet`](https://glmnet.stanford.edu/reference/glmnet.html) or
  [`ncvreg`](https://pbreheny.github.io/ncvreg/reference/ncvreg.html).

- `scaling`:

  The indices of the scale parameters. See `link{prestogp_fit}`.

- `nscale`:

  The number of scale parameters in the model.

- `common_scale`:

  Do all columns of locs have the same scale parameter?

- `param_sequence`:

  Records the indices of the various Matern parameters. See
  [`create_param_sequence`](https://niehs.github.io/PrestoGP/reference/create_param_sequence.md).

## References

- Apanasovich, T.V., Genton, M.G. and Sun, Y. "A valid Matérn class of
  cross-covariance functions for multivariate random fields with any
  number of components", Journal of the American Statistical
  Association (2012) 107(497):180-193.

- Messier, K.P. and Katzfuss, M. "Scalable penalized spatiotemporal
  land-use regression for ground-level nitrogen dioxide", The Annals of
  Applied Statistics (2021) 15(2):688-710.

## See also

[`VecchiaModel-class`](https://niehs.github.io/PrestoGP/reference/VecchiaModel-class.md),
[`FullModel-class`](https://niehs.github.io/PrestoGP/reference/FullModel-class.md),
[`MultivariateVecchiaModel-class`](https://niehs.github.io/PrestoGP/reference/MultivariateVecchiaModel-class.md),
[`prestogp_fit`](https://niehs.github.io/PrestoGP/reference/prestogp_fit-PrestoGPModel-method.md)

## Examples

``` r
pgp.vmodel <- new("VecchiaModel", n_neighbors = 25)
pgp.fmodel <- new("FullModel")
pgp.mmodel <- new("MultivariateVecchiaModel", n_neighbors = 25)
```
