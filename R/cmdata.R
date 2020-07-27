#' Simulated Crosswise Survey Data
#'
#' Small, artificially generated toy data set that comes in a cross-sectional
#' format where the unit of analysis is either country-year or
#' country-year-month. It provides artificial information for five countries
#' (Angola, Benin, France, Rwanda, and the UK) for a time span from 1990 to 1999 to
#' illustrate the use of the package.
#'
#' @docType data
#'
#' @usage data(cmdata)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{Y}{A binary indicator for the crosswise item (TRUE-TRUE or FALSE-FALSE) in the crosswise question}
#'  \item{A}{A binary indicator for the crosswise item (TRUE-TRUE or FALSE-FALSE) in the anchor question}
#'  \item{p}{An auxiliary probability in the crosswise question}
#'  \item{p.prime}{An auxiliary probability in the anchor question}
#' }
#' @references This data set was artificially created for the overviewR package.
#' @keywords datasets
#' @examples
#'
#' data(cmdata)
#' head(cmdata)
#'
"cmdata"
