#' Soil data from North Bavaria, Germany
#'
#' Soil physical and chemical data collected on a field in the Weissenstaedter
#' Becken, Germany
#'
#' @usage data(soil)
#' @format A data frame with the following columns:
#' \describe{
#' \item{x.coord:}{x coordinate (cm)}
#' \item{y.coord:}{y coordinate (cm)}
#' \item{nr:}{Sample number. The order of the sample number corresponds to the
#' order that the samples were collected.}
#' \item{moisture:}{Moisture content (Kg/Kg * 100\%)}
#' \item{NO3.N:}{Nitrate nitrogen (mg/Kg)}
#' \item{Total.N:}{Total nitrogen (mg/Kg)}
#' \item{NH4.N:}{Ammonium nitrogen (mg/Kg)}
#' \item{DOC:}{Dissolved organic carbon (mg/Kg)}
#' \item{N20N:}{Nitrous oxide (mg/Kg dried substance)}
#' }
#'
#' @docType data
#' @details For technical reasons some of the data were obtained as
#' differences of two measurements (which are not available anymore).
#' Therefore, some of the data have negative values.
#'
#' @author Martin Schlather, School of Business Informatics and Mathematics,
#' University of Mannheim
#'
#' @source The data were collected by Wolfgang Falk, Soil Physics Group,
#' University of Bayreuth, Germany. This data set was previously published in
#' the "RandomFields" R package.
#'
#' @references
#' \itemize{
#' \item Falk, W. "Kleinskalige raeumliche Variabilitaet von Lachgas und
#' bodenchemischen Parameters (Small scale spatial variability of nitrous
#' oxide and pedo-chemical parameters)" (2000). Master thesis, University of
#' Bayreuth, Germany.
#' }
"soil"

#' Pressure and temperature forecast errors over the Pacific Northwest
#'
#' This is a meteorological dataset that consists of differences between
#' forecasted and observed values of temperature and pressure at 157 locations
#' in the North American Pacific Northwest.
#'
#' @usage data(weather)
#'
#' @format A data frame with the following columns:
#' \describe{
#' \item{pressure:}{Difference between forecasted and observed pressure (Pa)}
#' \item{temperature:}{Difference between forecasted and observed temperature
#' (degrees C)}
#' \item{lon:}{Longitude coordinate of the location}
#' \item{lat:}{Latitude coordinate of the location}
#' }
#'
#' @docType data
#' @details These forecasts are from the University of Washington regional
#' numerical weather prediction ensemble (UWME; Grimit and Mass 2002;
#' Eckel and Mass 2005); they were valid on December 18, 2003 at 4 pm local
#' time, with a forecast horizon of 48 hours.
#'
#' @author Martin Schlather, School of Business Informatics and Mathematics,
#' University of Mannheim
#'
#' @source The data were obtained from Cliff Mass and Jeff Baars from the
#' University of Washington Department of Atmospheric Sciences. This data set
#' was previously published in the "RandomFields" R package.
#'
#' @references
#' \itemize{
#' \item Eckel, A. F. and Mass, C. F. "Aspects of effective mesoscale,
#' short-range ensemble forecasting", Weather and Forecasting (2005)
#' 20(3):328-350.
#' \item Gneiting, T., Kleiber, W. and Schlather, M. "Matern cross-covariance
#' functions for multivariate random fields", Journal of the American
#' Statistical Association (2010) 105(491):1167-1177.
#' \item Grimit, E. P. and Mass, C. F. "Initial results of a mesoscale
#' short-range forecasting system over the Pacific Northwest", Weather
#' and Forecasting 17(2):192-205.
#' }
"weather"
