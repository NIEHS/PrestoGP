# PrestoGP 0.2.0.9039 (2024-12-26)

## MINOR IMPROVEMENTS

  * Rewrote some tests that were failing on some verions of R

# PrestoGP 0.2.0.9038 (2024-12-23)

## NEW FEATURES

  * Added the method `get_Y` to extract the Y values (including the imputed
    Y values) from a PrestoGP model

  * Added the method `get_linear_model` to extract the fitted cv.glmnet object
    from a PrestoGP model
	
  * Added the method `plot_beta` to plot the glide path of the
    regression coefficients

## MINOR IMPROVEMENTS

  * Added a NEWS.md file to track changes in each update

  * Modified the `show` method to not print zero betas
  
  * Modified `prestogp_predict` to retun predictions in the form of a
    list for multivariate models

## BUG FIXES

  * Fixed a bug in the intercept calculation in `get_beta`

  * Fixed a bug where multivariate models with no variable names for X were
    being assigned arbitrary variable names

## DOCUMENTATION FIXES

  * Updated the documention in `PrestoGPModel-class` to reflect that the
    Y_train slot is updated to contain imputed values during model
    fitting
	
  * Added documentation for all new methods
