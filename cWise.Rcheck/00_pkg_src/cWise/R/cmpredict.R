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
#' @param nsim Number of parametric-bootstrap draws.
#' @param seed Optional integer seed for reproducible bootstrap draws. When set,
#' the caller's RNG state is restored before returning.
#' @param draws If `TRUE`, attach the raw bootstrap draw matrix as a `"draws"`
#' attribute on the returned data frame.
#'
#' @return A data frame with `estimate`, `conf.low`, and `conf.high` columns,
#' with one row for each prediction scenario. When `draws = TRUE`, its `"draws"`
#' attribute contains the raw parametric-bootstrap matrix.
#' @seealso [cmreg()] to fit the required outcome-model object and
#' [cmpredict_p()] for predictions from a predictor model.
#' @references Atsusaka, Y. and Stevenson, R. T. (2023). The crosswise model
#' for sensitive survey questions. \doi{10.1017/pan.2021.43}.
#' @examples
#' m <- cmreg(Y ~ female + age, anchor = A, p = 0.1, p.prime = 0.15,
#'            data = cmdata2)
#' predictions <- cmpredict(m, typical = c(age = 30),
#'                           zval = c(female = 0, female = 1))
#' predictions
#' @export
#' @importFrom mvtnorm "rmvnorm"


cmpredict <- function(out, newdata = NULL, zval = NULL, typical = NULL,
                      nsim = 1000L, seed = NULL, draws = FALSE){

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
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L || nsim != as.integer(nsim)) {
    stop("`nsim` must be a single positive integer.", call. = FALSE)
  }
  if (!is.logical(draws) || length(draws) != 1L || is.na(draws)) {
    stop("`draws` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

# GRAB FULL-PRECISION COEFFICIENTS
  idx <- cm_par_index(ncol(typ.vec), model = "outcome")
  coef.beta  = out$estimates[idx$beta]
  vcovs = out$VCV[idx$beta, idx$beta, drop = FALSE]

# PARAMETRIC BOOTSTRAP
  coef.sim <- cm_with_seed(seed, function() {
    mvtnorm::rmvnorm(n = as.integer(nsim), mean = coef.beta, sigma = vcovs)
  })

  lin.agg <- typ.vec %*% t(coef.sim) # Linear aggregator
  pi.sim = i.logit(lin.agg) # Inverse logit

  point_estimate <- i.logit(typ.vec %*% coef.beta)
  cm_prediction_summary(pi.sim, point_estimate, draws)
}
