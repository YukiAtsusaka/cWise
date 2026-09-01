#' Print and summarise crosswise regression fits
#'
#' @param x A `cmreg` or `cmreg_p` fit object.
#' @param object A `cmreg` or `cmreg_p` fit object.
#' @param digits Number of significant digits to display.
#' @param ... Additional arguments, currently ignored.
#'
#' @return `summary.cmreg()` returns a list containing the fit call, coefficient
#' tables, log-likelihood, sample size, convergence code, and model type.
#'
#' @name cmreg-methods
NULL

#' @rdname cmreg-methods
#' @export
print.cmreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n")
  print(x$Call)
  cat("\nCrosswise regression with the latent trait as ",
      if (identical(x$model, "predictor")) "a predictor" else "an outcome", ".\n",
      sep = "")
  cat("n = ", x$n, ", log-likelihood = ", format(x$logLik, digits = digits), "\n\n", sep = "")

  print(summary(x), digits = digits)
  invisible(x)
}

#' @rdname cmreg-methods
#' @export
summary.cmreg <- function(object, ...) {
  result <- list(
    call = object$Call,
    coefficients = object$Coefficients,
    auxiliary = object$AuxiliaryCoef,
    logLik = object$logLik,
    n = object$n,
    convergence = object$convergence,
    model = object$model
  )
  if (inherits(object, "cmreg_p")) {
    result$auxiliary2 <- object$AuxiliaryCoef2
  }
  class(result) <- "summary.cmreg"
  result
}

#' @rdname cmreg-methods
#' @export
print.summary.cmreg <- function(x,
                                digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("Coefficients:\n")
  print(round(x$coefficients, digits = digits))
  cat("\nAuxiliary coefficients:\n")
  print(round(x$auxiliary, digits = digits))
  if (!is.null(x$auxiliary2)) {
    cat("\nAdditional auxiliary coefficients:\n")
    print(round(x$auxiliary2, digits = digits))
  }
  if (!identical(x$convergence, 0L)) {
    cat("\nOptimization convergence code: ", x$convergence, "\n", sep = "")
  }
  invisible(x)
}
