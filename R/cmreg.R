#' @title cmreg
#'
#' @description \code{cmreg} is used to run a regression with the latent sensitive trait as an outcome.
#'
#' @param formula an object of class "formula":
#' a symbolic description of the model to be fitted.
#' Ex. Crosswise response ~ Covariates. The anchor response is supplied
#' separately through `anchor`.
#' @param anchor Unquoted column name (or a single character column name) for the
#' anchor response.
#' @param p an auxiliary probability for the crosswise question.
#' @param p.prime an auxiliary probability for the anchor question.
#' @param data a data frame containing information from the crosswise model and covariates.
#' @param start Optional numeric vector of starting values for the beta and theta
#' parameters. By default, starts are derived from binomial GLMs for the observed
#' crosswise and anchor responses.
#' @param n.start Number of optimization starts, including the data-informed start.
#' @param control A list of control settings passed to [stats::optim()]. `fnscale`
#' is fixed internally because the log-likelihood is maximized.
#'
#' @return An object of class `cmreg`, a list with:
#' \describe{
#'   \item{Call}{The matched model call.}
#'   \item{Coefficients}{A coefficient table for the latent-trait outcome model.}
#'   \item{AuxiliaryCoef}{A coefficient table for the anchor-response model.}
#'   \item{VCV}{The estimated variance-covariance matrix of all parameters.}
#'   \item{estimates, std.errors}{Named full-precision parameter estimates and
#'   standard errors.}
#'   \item{logLik, n, convergence}{The optimized log-likelihood, number of
#'   complete observations, and optimizer convergence code.}
#' }
#' @seealso [cmpredict()] for predicted latent-trait probabilities and
#' [cmreg_p()] for a model with the latent trait as a predictor.
#' @references Atsusaka, Y. and Stevenson, R. T. (2023). The crosswise model
#' for sensitive survey questions. \doi{10.1017/pan.2021.43}.
#' @examples
#' m <- cmreg(Y ~ female + age, anchor = A, p = 0.1, p.prime = 0.15,
#'            data = cmdata2)
#' m
#' @export


cmreg <- function(formula, anchor, p, p.prime, data, start = NULL,
                  n.start = 3L, control = list()){

  call <- match.call()

  anchor_info <- cm_data_column(data, substitute(anchor), "anchor", parent.frame())
  ## `model.frame()` evaluates transformations such as poly() before applying
  ## its `na.action`.  Remove rows missing raw formula variables first so that
  ## transformed formulas handle incomplete data consistently.
  formula_complete <- complete.cases(data[, all.vars(formula), drop = FALSE])
  mf <- model.frame(formula, data[formula_complete, , drop = FALSE],
                    na.action = na.fail)
  anchor_values <- anchor_info$value[formula_complete]
  keep <- complete.cases(anchor_values)
  df <- mf[keep, , drop = FALSE]
  X1 <- model.matrix(formula, df)
  A <- anchor_values[keep]
  Y <- model.response(df)
  k <- dim(X1)[2]      # Number of beta parameters
  idx <- cm_par_index(k, model = "outcome")

  cm_validate_probabilities(p, p.prime)
  cm_validate_binary(Y, "formula response")
  cm_validate_binary(A, "anchor")
  if (nrow(X1) <= idx$npar) {
    stop("The number of complete observations must exceed the number of parameters.",
         call. = FALSE)
  }
  cm_validate_no_separation(X1, Y, "formula response")
  cm_validate_no_separation(X1, A, "anchor")

  init <- if (is.null(start)) {
    c(cm_glm_start(X1, Y), cm_glm_start(X1, A))
  } else {
    as.numeric(start)
  }
  if (length(init) != idx$npar || any(!is.finite(init))) {
    stop(sprintf("`start` must contain %d finite beta/theta values.", idx$npar),
         call. = FALSE)
  }

  # LOG-LIKELIHOOD FUNCTION
  log.L <- function(par) {
    cm_outcome_loglik(par, X1, Y, A, p, p.prime)
  }


  # MAXIMIZATION
  MLE = cm_multistart_optim(log.L, init, n.start, control)
  H = cm_outcome_hessian(MLE$par, X1, Y, A, p, p.prime)
  VCV = cm_hessian_vcov(H, idx$npar)
  SE = sqrt(diag(VCV))


# OUTPUT

  Mlist <- list()

  z = MLE$par / SE
  pv = 2*(1- pnorm(abs(z)))
  Mlist[[1]] <- call
  Mlist[[2]] <- t(rbind(MLE$par[idx$beta], SE[idx$beta], z[idx$beta], pv[idx$beta]))
  Mlist[[3]] <- t(rbind(MLE$par[idx$theta], SE[idx$theta],
                        z[idx$theta], pv[idx$theta]))
  Mlist[[4]] <- VCV # Estimated variance-covariance matrix
  colnames(Mlist[[2]]) <- c("Estimate", "Std. Error", "z score", "Pr(>|z|)")
  colnames(Mlist[[3]]) <- c("Estimate", "Std. Error", "z score", "Pr(>|z|)")

  varnam <- colnames(X1)
  rownames(Mlist[[2]]) <- varnam
  rownames(Mlist[[3]]) <- varnam

  parameter_names <- c(paste0("beta:", varnam), paste0("theta:", varnam))
  names(MLE$par) <- parameter_names
  names(SE) <- parameter_names
  dimnames(Mlist[[4]]) <- list(parameter_names, parameter_names)

  names(Mlist) <- c("Call", "Coefficients", "AuxiliaryCoef", "VCV")
  Mlist$estimates <- MLE$par
  Mlist$std.errors <- SE
  Mlist$logLik <- MLE$value
  Mlist$hessian <- H
  Mlist$n <- length(Y)
  Mlist$convergence <- MLE$convergence
  Mlist$p <- p
  Mlist$p.prime <- p.prime
  Mlist$model <- "outcome"
  # The model-frame terms retain transformation metadata (for example, the
  # orthogonal-polynomial coefficients needed to predict from `poly()`).
  Mlist$terms <- terms(mf)
  Mlist$xlevels <- stats::.getXlevels(Mlist$terms, df)
  Mlist["contrasts"] <- list(attr(X1, "contrasts"))
  Mlist$design_columns <- varnam
  class(Mlist) <- "cmreg"


  return(Mlist)

}
