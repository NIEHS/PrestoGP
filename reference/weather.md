# Pressure and temperature forecast errors over the Pacific Northwest

This is a meteorological dataset that consists of differences between
forecasted and observed values of temperature and pressure at 157
locations in the North American Pacific Northwest.

## Usage

``` r
data(weather)
```

## Format

A data frame with the following columns:

- pressure::

  Difference between forecasted and observed pressure (Pa)

- temperature::

  Difference between forecasted and observed temperature (degrees C)

- lon::

  Longitude coordinate of the location

- lat::

  Latitude coordinate of the location

## Source

The data were obtained from Cliff Mass and Jeff Baars from the
University of Washington Department of Atmospheric Sciences. This data
set was previously published in the "RandomFields" R package.

## Details

These forecasts are from the University of Washington regional numerical
weather prediction ensemble (UWME; Grimit and Mass 2002; Eckel and Mass
2005); they were valid on December 18, 2003 at 4 pm local time, with a
forecast horizon of 48 hours.

## References

- Eckel, A. F. and Mass, C. F. "Aspects of effective mesoscale,
  short-range ensemble forecasting", Weather and Forecasting (2005)
  20(3):328-350.

- Gneiting, T., Kleiber, W. and Schlather, M. "Matern cross-covariance
  functions for multivariate random fields", Journal of the American
  Statistical Association (2010) 105(491):1167-1177.

- Grimit, E. P. and Mass, C. F. "Initial results of a mesoscale
  short-range forecasting system over the Pacific Northwest", Weather
  and Forecasting 17(2):192-205.

## Author

Martin Schlather, School of Business Informatics and Mathematics,
University of Mannheim
