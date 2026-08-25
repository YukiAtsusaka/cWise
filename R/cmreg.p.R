#' @title cmreg.p
#'
#' @description \code{cmreg} is used to run a regression with the latent sensitive trait as a predictor.
#'
#' @param formula an object of class "formula":
#' a symbolic description of the model to be fitted.
#' Ex. Outcome ~ Covariates. The crosswise and anchor responses are supplied
#' separately through `crosswise` and `anchor`.
#' @param crosswise Unquoted column name (or a single character column name) for
#' the crosswise response.
#' @param anchor Unquoted column name (or a single character column name) for the
#' anchor response.
#' @param p an auxiliary probability for the crosswise question.
#' @param p.prime an auxiliary probability for the anchor question.
#' @param data a data frame containing information from the crosswise model, the outcome variable, and covariates.
#' @param start Optional numeric vector of starting values on the reporting scale:
#' beta, theta, gamma, crosswise effect, and a positive sigma. By default, starts
#' are derived from GLM and least-squares fits to observed responses.
#' @param n.start Number of optimization starts, including the data-informed start.
#' @param control A list of control settings passed to [stats::optim()]. `fnscale`
#' is fixed internally because the log-likelihood is maximized.
#'
#' @return A list containing the estimated results and related statistics.
#' @examples
#' m2 <- cmreg.p(V ~ age + female, crosswise = Y, anchor = A, p = 0.1,
#'               p.prime = 0.15, data = cmdata3)
#' m2
#' @export


cmreg.p <- function(formula, crosswise, anchor, p, p.prime, data, start = NULL,
                    n.start = 3L, control = list()){
  call <- match.call()
  i.logit <- function(XB){ exp(XB)/(1 + exp(XB))}
  crosswise_info <- cm_data_column(data, substitute(crosswise), "crosswise", parent.frame())
  anchor_info <- cm_data_column(data, substitute(anchor), "anchor", parent.frame())
  mf <- model.frame(formula, data, na.action = na.pass)
  keep <- complete.cases(mf, crosswise_info$value, anchor_info$value)
  df <- mf[keep, , drop = FALSE]
  X1 <- model.matrix(formula, df)

  A <- anchor_info$value[keep]
  Y <- crosswise_info$value[keep]
  V <- model.response(df)

  k <- dim(X1)[2]      # Number of beta parameters ( #covariate + 1)
  idx <- cm_par_index(k, model = "predictor")

  init <- if (is.null(start)) {
    gamma <- cm_lm_start(X1, V)
    residual_sd <- sqrt(mean((V - drop(X1 %*% gamma))^2))
    if (!is.finite(residual_sd) || residual_sd <= 0) residual_sd <- 1
    c(cm_glm_start(X1, Y), cm_glm_start(X1, A), gamma, 0, residual_sd)
  } else {
    as.numeric(start)
  }
  if (length(init) != idx$npar || any(!is.finite(init))) {
    stop(sprintf("`start` must contain %d finite parameter values.", idx$npar),
         call. = FALSE)
  }
  if (init[idx$sigma] <= 0) {
    stop("The sigma element of `start` must be positive.", call. = FALSE)
  }
  init[idx$log_sigma] <- log(init[idx$sigma])

#---------------------------------------------------------------------------------------#
# PARAMTERS TO ESTIMATE (11 parameters, if two covariates)
# beta0, beta1, beta2,
# theta0, theta1, theta2,
# gamma0, gamma1, gamma2, gamma.cm (Z), sigma
#---------------------------------------------------------------------------------------#

# LOG-LIKELIHOOD FUNCTION
log.L.pred <- function(par) {
  sigma <- exp(par[idx$log_sigma])

  sum(log( (1/(sigma*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[idx$gamma] + 1*par[idx$gamma_z] ))^2/(2*(sigma^2)) )
           * (i.logit(X1 %*% par[idx$beta])) * (i.logit(X1 %*% par[idx$theta]))  * (p^Y)*((1-p)^(1-Y)) * ((1-p.prime)^A)*(p.prime^(1-A))

        +  (1/(sigma*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[idx$gamma] + 0*par[idx$gamma_z] ))^2/(2*(sigma^2)) )
           * (1-(i.logit(X1 %*% par[idx$beta]))) * (i.logit(X1 %*% par[idx$theta])) * ((1-p)^Y)*(p^(1-Y))*((1-p.prime)^A)*(p.prime^(1-A))

        +  (1/(sigma*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[idx$gamma] + 1*par[idx$gamma_z] ))^2/(2*(sigma^2)) )
           * (i.logit(X1 %*% par[idx$beta])) * (1 - i.logit(X1 %*% par[idx$theta])) * (1/4)

        +  (1/(sigma*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[idx$gamma] + 0*par[idx$gamma_z] ))^2/(2*(sigma^2)) )
           * (1-(i.logit(X1 %*% par[idx$beta]))) * (1 - i.logit(X1 %*% par[idx$theta])) * (1/4)
  ))
}


# MAXIMIZATION
MLE = cm_multistart_optim(log.L.pred, init, n.start, control)
VCV = cm_hessian_vcov(MLE$hessian, idx$npar) # VCV on the optimization scale
estimate = MLE$par
estimate[idx$sigma] = exp(MLE$par[idx$log_sigma])

# Delta-transform the covariance matrix from log(sigma) back to sigma for
# reporting. This leaves all other parameter variances and covariances intact.
jacobian = diag(idx$npar)
jacobian[idx$sigma, idx$log_sigma] = estimate[idx$sigma]
VCV = jacobian %*% VCV %*% jacobian
SE = sqrt(diag(VCV))


# OUTPUT

  Mlist <- list()
  z = estimate / SE
  pv = 2*(1- pnorm(abs(z)))
  Mlist[[1]] <- call
  gamma_idx <- c(idx$gamma, idx$gamma_z)
  theta_sigma_idx <- c(idx$theta, idx$sigma)
  Mlist[[2]] <- t(rbind(estimate[gamma_idx],
                             SE[gamma_idx],
                              z[gamma_idx],
                             pv[gamma_idx]))
  Mlist[[3]] <- t(rbind(estimate[idx$beta], SE[idx$beta],
                               z[idx$beta], pv[idx$beta]))
  Mlist[[4]] <- t(rbind(estimate[theta_sigma_idx],
                             SE[theta_sigma_idx],
                              z[theta_sigma_idx],
                             pv[theta_sigma_idx]))
  colnames(Mlist[[2]]) <- c("Estimate", "Std. Error", "Z score", "Pr(>|z|)")
  colnames(Mlist[[3]]) <- c("Estimate", "Std. Error", "Z score", "Pr(>|z|)")
  colnames(Mlist[[4]]) <- c("Estimate", "Std. Error", "Z score", "Pr(>|z|)")
  Mlist[[5]] <- VCV # Estimated variance-covariance matrix on the reporting scale

  varnam <- colnames(X1)
  rownames(Mlist[[2]]) <- c(varnam, crosswise_info$name)
  rownames(Mlist[[3]]) <- varnam
  rownames(Mlist[[4]]) <- c(varnam, "sigma")

  parameter_names <- c(
    paste0("beta:", varnam),
    paste0("theta:", varnam),
    paste0("gamma:", varnam),
    paste0("gamma_z:", crosswise_info$name),
    "sigma"
  )
  names(estimate) <- parameter_names
  names(SE) <- parameter_names
  dimnames(Mlist[[5]]) <- list(parameter_names, parameter_names)

names(Mlist) <- c("Call", "Coefficients", "AuxiliaryCoef", "AuxiliaryCoef2", "VCV")
Mlist$estimates <- estimate
Mlist$std.errors <- SE
Mlist$logLik <- MLE$value
Mlist$n <- length(V)
Mlist$convergence <- MLE$convergence
Mlist$p <- p
Mlist$p.prime <- p.prime
Mlist$model <- "predictor"
Mlist$terms <- terms(mf)
Mlist$xlevels <- stats::.getXlevels(Mlist$terms, df)
Mlist["contrasts"] <- list(attr(X1, "contrasts"))
Mlist$design_columns <- varnam
class(Mlist) <- c("cmreg.p", "cmreg")

return(Mlist)
}
