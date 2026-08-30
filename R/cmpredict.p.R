#' @title cmpredict_p
#'
#' @description Perform a post-estimation prediction with uncertainty quantification via parametric bootstrap
#'
#' @param out An output of \code{cmreg_p}.
#' @param newdata A named data frame, list, or vector supplying values for every
#' covariate in the fitted formula. A data frame returns absent/present latent-
#' trait scenarios for every row.
#' @param zval Optional named vector for one covariate to vary. Its name must
#' match a formula variable; all remaining covariates are supplied through
#' `newdata` or `typical`.
#' @param typical Optional named vector or list of fixed covariate values. This
#' is a convenience alternative to `newdata` when used with `zval`.
#' @param type Prediction scale: `"response"` (default) or `"link"`. The
#' Gaussian outcome model has an identity link, so these are currently equal.
#' @param nsim Number of parametric-bootstrap draws.
#' @param seed Optional integer seed for reproducible bootstrap draws. When set,
#' the caller's RNG state is restored before returning.
#' @param draws If `TRUE`, attach the raw bootstrap draw matrix as a `"draws"`
#' attribute on the returned data frame.
#'
#' @return A data frame with `estimate`, `conf.low`, and `conf.high` columns.
#' Each input scenario contributes an absent-trait row followed by a present-
#' trait row. When `draws = TRUE`, its `"draws"` attribute contains the raw
#' parametric-bootstrap matrix.
#' @examples
#' m2 <- cmreg_p(V ~ age + female, crosswise = Y, anchor = A, p = 0.1,
#'               p.prime = 0.15, data = cmdata3)
#' predictions <- cmpredict_p(m2, newdata = data.frame(age = 30, female = 1))
#' predictions
#' @export
#' @importFrom mvtnorm "rmvnorm"

cmpredict_p <- function(out, newdata = NULL, zval = NULL, typical = NULL,
                        type = c("response", "link"), nsim = 1000L,
                        seed = NULL, draws = FALSE){

  type <- match.arg(type)

# VALIDATE FIT AND BUILD DESIGN MATRIX
  if (!inherits(out, "cmreg_p") || !identical(out$model, "predictor")) {
    stop("`out` must be a fitted predictor-model `cmreg_p` object.", call. = FALSE)
  }
  required_fields <- c("terms", "design_columns", "estimates")
  if (!all(required_fields %in% names(out))) {
    stop("`out` was created by an older cWise version; refit the model before predicting.",
         call. = FALSE)
  }
  covariate_design <- cm_prediction_data(out, newdata, zval, typical)
  if (length(nsim) != 1L || is.na(nsim) || nsim < 1L || nsim != as.integer(nsim)) {
    stop("`nsim` must be a single positive integer.", call. = FALSE)
  }
  if (!is.logical(draws) || length(draws) != 1L || is.na(draws)) {
    stop("`draws` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

# GRAB FULL-PRECISION COEFFICIENTS AND MATCHED COVARIANCE MATRIX
  idx <- cm_par_index(ncol(covariate_design), model = "predictor")
  gamma_idx <- c(idx$gamma, idx$gamma_z)
  coef.gamma <- out$estimates[gamma_idx]
  vcovs = out$VCV[gamma_idx, gamma_idx, drop = FALSE]


# PARAMETRIC BOOTSTRAP
  coef.sim <- cm_with_seed(seed, function() {
    rmvnorm(n = as.integer(nsim), mean = coef.gamma, sigma = vcovs)
  })


# ABSENT/PRESENT LATENT-TRAIT SCENARIOS
  scenario_index <- rep(seq_len(nrow(covariate_design)), each = 2L)
  latent_trait <- rep(c(0, 1), times = nrow(covariate_design))
  typ.vec <- cbind(covariate_design[scenario_index, , drop = FALSE], latent_trait)
  rownames(typ.vec) <- paste0(
    rownames(covariate_design)[scenario_index],
    ifelse(latent_trait == 0, ": trait absent", ": trait present")
  )

  lin.agg <- typ.vec %*% t(coef.sim) # Linear aggregator

  # The current Gaussian model uses an identity link. Keeping the switch here
  # establishes the same contract as predict.glm() for future link functions.
  point_estimate <- typ.vec %*% coef.gamma
  prediction <- switch(type, response = lin.agg, link = lin.agg)
  cm_prediction_summary(prediction, point_estimate, draws)
}

#' @rdname cmpredict_p
#' @param ... Arguments passed to \code{cmpredict_p()}.
#' @export
cmpredict.p <- function(...) {
  .Deprecated("cmpredict_p")
  cmpredict_p(...)
}
