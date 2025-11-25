# PrestoGP 0.2.0.9050 (2025-11-13)

## BUG FIXES

* Fixed a bug where covariance matrices were computed incorrectly during
  imputation

* Fixed a bug that caused the imputation procedure to crash for univariate
  models without a common scale parameter

* Fixed a bug where the incorrect number of imputations was performed

## MINOR IMPROVEMENTS

* Added additional tests to confirm that imputation covariance matrices are
  computed correctly

# PrestoGP 0.2.0.9049 (2025-5-13)

## BUG FIXES

* Fixed another bug in `prestogp_predict` caused by numerically singular
  matrices

# PrestoGP 0.2.0.9048 (2025-5-2)

## BUG FIXES

* Fixed a bug that occurred when lod.upper or lod.lower is a vector of values

## MINOR IMPROVEMENTS

* Modified some tests to verify that the bug mentioned above is patched
  and does not recur

# PrestoGP 0.2.0.9047 (2025-4-28)

## BUG FIXES

* Fixed a bug in the multivariate version of `prestogp_predict`

# PrestoGP 0.2.0.9046 (2025-4-24)

## BREAKING CHANGES

* The parameter `lod` has been renamed to be `lod.upper`

* The parameter `max_iters` has been renamed as `max.iters` for a consistent
  naming convention

## NEW FEATURES

* The `prestogp_fit` method now supports lower limits of detection via the
  `lod.lower` parameter

* Additional tuning parameters can be specified for the LOD imputation
  procedure, namely `n.impute`, `eps.impute`, and `maxit.impute`

## BUG FIXES

* Added another set of fixes to avoid crashes when the estimated U matrix is
  numerically singular

* Subsampling is once again used for initial scale parameter estimates for
  large (>10,000 observation) data sets for computational reasons

# PrestoGP 0.2.0.9045 (2025-3-10)

## MINOR IMPROVEMENTS

* Fixed some linting issues

# PrestoGP 0.2.0.9044 (2025-3-7)

## BREAKING CHANGES

* Removed the `relax` parameter from `prestogp_fit`

* All glmnet-based models now use lambda.min (rather than lambda.1se) to
  choose the optimal value of lambda

## NEW FEATURES

* Added the parameter `penalty` to `prestogp_fit` to specify the type of
  penalized regression (relaxed lasso, SCAD, and MCP)

* Added the parameter `alpha` to `prestogp_fit` to change the alpha parameter
  in glmnet/ncvreg (which can be used to fit ridge regression or elastic
  net models)

## MINOR IMPROVEMENTS

* Subsampling is no longer used to calculate the initial estimates of the
  scale parameters in `calc_covparams`

* Fixed an issue where the C++ code displayed repeated warnings when an
  estimated submatrix was nearly numerically singular

* Many new tests were added to test the new penalties

## BUG FIXES

* Fixed an issue where initial scale parameter estimates were sometimes set
  to 0 due to subsampling

* Fixed an issue where likelihood calculations could crash if the estimated
  U matrix is numerically singular

* Fixed an issue where imputed values were not being updated properly

# PrestoGP 0.2.0.9043 (2025-1-29)

## BUG FIXES

* Imputation code has been rewritten to avoid computing the full
  covariance matrix, resulting in a huge performance improvement
  
* The nntmvn package dependency has been removed to fix errors caused
  by a recent update

# PrestoGP 0.2.0.9042 (2025-1-14)

## BREAKING CHANGES

* The slot `lambda_1se_idx` has been removed from the `PrestoGPModel`
  superclass

## NEW FEATURES

* Added the option `relax` to `prestogp_fit` to allow fitting a
  relaxed lasso model

## MINOR IMPROVEMENTS

* Eliminated the type checking on the `linear_model` slot in the
  `PrestoGPModel` superclass, which will allow new types of regression
  models (e.g., relaxed lasso, SCAD, MCP)

# PrestoGP 0.2.0.9041 (2025-1-7)

## NEW FEATURES

* Added the option `adaptive` to `prestogp_fit` to allow fitting an
  adaptive lasso model

## MINOR IMPROVEMENTS

* Eliminated some redundant code in `estimate_betas`

# PrestoGP 0.2.0.9040 (2025-1-3)

## BREAKING CHANGES

* The function `create.param.sequence` has been renamed to be
  `create_param_sequence` for consistency with other function names

* The parameter `apanasovich` has been renamed `common_scale`
  throughout the package
  
## MINOR IMPROVEMENTS

* Tests were updated to reflect the new `common_scale` parameter and
  the `create_param_sequence` function

* Fixed some formatting issues in NEWS.md

## DOCUMENTATION FIXES

* The documentation for `create_param_seqence` was updated to clarify
  how to extract/assign scale parameters when there is more than one
  scale parameter per outcome

* The documentation for `PrestoGP-Model-class` and `prestogp_fit` has
  been updated to reflect the new `common_scale` parameter

# PrestoGP 0.2.0.9039 (2024-12-26)

## MINOR IMPROVEMENTS

* Rewrote some tests that were failing on some verions of R

# PrestoGP 0.2.0.9038 (2024-12-23)

## NEW FEATURES

* Added the method `get_Y` to extract the Y values (including the
  imputed Y values) from a PrestoGP model 

* Added the method `get_linear_model` to extract the fitted cv.glmnet
  object from a PrestoGP model
	
* Added the method `plot_beta` to plot the glide path of the
  regression coefficients 

## MINOR IMPROVEMENTS

* Added a NEWS.md file to track changes in each update

* Modified the `show` method to not print zero betas
  
* Modified `prestogp_predict` to retun predictions in the form of a
  list for multivariate models

## BUG FIXES

* Fixed a bug in the intercept calculation in `get_beta`

* Fixed a bug where multivariate models with no variable names for X
  were being assigned arbitrary variable names

## DOCUMENTATION FIXES

* Updated the documention in `PrestoGPModel-class` to reflect that the
  Y_train slot is updated to contain imputed values during model fitting
	
* Added documentation for all new methods
