#' Simulated Crosswise Survey Data III
#'
#' Typical, artificially generated toy data set that comes with the binary response in the crosswise question (Y=1 if TRUE-TRUE or FALSE-FALSE, Y=0 otherwise),
#' the binary response in the anchor question, the outcome variable, two covariates, and auxiliary probabilities in the crosswise and anchor questions, respectively (these columns are only included here for an illustrative purpose).
#' @docType data
#'
#' @usage data(cmdata3)
#'
#' @format an object of class \code{"data.frame"}
#' \describe{
#'  \item{V}{an outcome variable}
#'  \item{female}{A binary covariate representing being a female}
#'  \item{age}{A discrete covariate representing age}
#'  \item{Y}{a binary indicator for the crosswise item (TRUE-TRUE or FALSE-FALSE) in the crosswise question}
#'  \item{A}{a binary indicator for the crosswise item (TRUE-TRUE or FALSE-FALSE) in the anchor question}
#'  \item{p}{an auxiliary probability in the crosswise question}
#'  \item{p.prime}{an auxiliary probability in the anchor question}
#' }
#' @references This data set was artificially created for the cWise package.
#' @keywords datasets
#' @examples
#'
#' data(cmdata3)
#' head(cmdata3)
#'
"cmdata3"
