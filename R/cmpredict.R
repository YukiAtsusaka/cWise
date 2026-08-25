#' @title cmpredict
#'
#' @description Perform a post-estimation prediction with uncertainty quantification via parametric bootstrap
#'
#' @param out An output of \code{cmreg}.
#' @param newdata A named data frame, list, or vector supplying values for every
#' covariate in the fitted formula. A data frame returns one prediction row per
#' row of `newdata`.
#' @param zval Optional named vector for one covariate to vary. Its name must
#' match a formula variable; all remaining covariates are supplied through
#' `newdata` or `typical`.
#' @param typical Optional named vector or list of fixed covariate values. This
#' is a convenience alternative to `newdata` when used with `zval`.
#'
#' @return A matrix of parametric-bootstrap draws of predicted probabilities,
#' with one row for each prediction scenario.
#' @examples
#' m <- cmreg(Y ~ female + age, anchor = A, p = 0.1, p.prime = 0.15,
#'            data = cmdata2)
#' predictions <- cmpredict(m, typical = c(age = 30),
#'                           zval = c(female = 0, female = 1))
#' predictions
#' @export
#' @importFrom mvtnorm "rmvnorm"


cmpredict <- function(out, newdata = NULL, zval = NULL, typical = NULL){

i.logit <- function(XB){ exp(XB)/(1 + exp(XB))}

# VALIDATE FIT AND BUILD DESIGN MATRIX
  if (!inherits(out, "cmreg") || identical(out$model, "predictor")) {
    stop("`out` must be a fitted outcome-model `cmreg` object.", call. = FALSE)
  }
  required_fields <- c("terms", "design_columns", "estimates")
  if (!all(required_fields %in% names(out))) {
    stop("`out` was created by an older cWise version; refit the model before predicting.",
         call. = FALSE)
  }
  typ.vec <- cm_prediction_data(out, newdata, zval, typical)

# GRAB FULL-PRECISION COEFFICIENTS
  idx <- cm_par_index(ncol(typ.vec), model = "outcome")
  coef.beta  = out$estimates[idx$beta]
  vcovs = out$VCV[idx$beta, idx$beta, drop = FALSE]

# PARAMETRIC BOOTSTRAP
  set.seed(20200730)
  coef.sim <- mvtnorm::rmvnorm(n=10000, mean=coef.beta, sigma=vcovs) # from mvtnorm

  lin.agg <- typ.vec %*% t(coef.sim) # Linear aggregator
  pi.sim = i.logit(lin.agg) # Inverse logit

  return(pi.sim)
}
