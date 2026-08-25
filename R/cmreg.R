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
#' @return A list containing the estimated results and related statistics.
#' @examples
#' m <- cmreg(Y ~ female + age, anchor = A, p = 0.1, p.prime = 0.15,
#'            data = cmdata2)
#' m
#' @export


cmreg <- function(formula, anchor, p, p.prime, data, start = NULL,
                  n.start = 3L, control = list()){

  call <- match.call()

  i.logit <- function(XB){ exp(XB)/(1 + exp(XB))}

  anchor_info <- cm_data_column(data, substitute(anchor), "anchor", parent.frame())
  mf <- model.frame(formula, data, na.action = na.pass)
  keep <- complete.cases(mf, anchor_info$value)
  df <- mf[keep, , drop = FALSE]
  X1 <- model.matrix(formula, df)
  A <- anchor_info$value[keep]
  Y <- model.response(df)
  k <- dim(X1)[2]      # Number of beta parameters
  idx <- cm_par_index(k, model = "outcome")

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
    sum(   Y*log(            ((2*p-1)*(i.logit(X1 %*% par[idx$beta])) + (0.5-p) )*i.logit(X1 %*% par[idx$theta]) + 0.5)
           + (1-Y)*log( 1 - (((2*p-1)*(i.logit(X1 %*% par[idx$beta])) + (0.5-p) )*i.logit(X1 %*% par[idx$theta]) + 0.5))
           + A*log(        (0.5-p.prime)*i.logit(X1 %*% par[idx$theta]) + 0.5 )
           + (1-A)*log(1- ((0.5-p.prime)*i.logit(X1 %*% par[idx$theta]) + 0.5 ))
    )
  }


  # MAXIMIZATION
  MLE = cm_multistart_optim(log.L, init, n.start, control)
  VCV = cm_hessian_vcov(MLE$hessian, idx$npar)
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
