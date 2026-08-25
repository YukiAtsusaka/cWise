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
#'
#' @return A list containing the estimated results and related statistics.
#' @examples
#' m2 <- cmreg.p(V ~ age + female, crosswise = Y, anchor = A, p = 0.1,
#'               p.prime = 0.15, data = cmdata3)
#' m2
#' @export


cmreg.p <- function(formula, crosswise, anchor, p, p.prime, data){
  options(warn=-1)

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

  init <- rep(0.01, idx$npar) # Initial values for optimization
  init[idx$sigma] <- 1

#---------------------------------------------------------------------------------------#
# PARAMTERS TO ESTIMATE (11 parameters, if two covariates)
# beta0, beta1, beta2,
# theta0, theta1, theta2,
# gamma0, gamma1, gamma2, gamma.cm (Z), sigma
#---------------------------------------------------------------------------------------#

# LOG-LIKELIHOOD FUNCTION
log.L.pred <- function(par) {
  sum(log( (1/(par[idx$sigma]*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[idx$gamma] + 1*par[idx$gamma_z] ))^2/(2*(par[idx$sigma]^2)) )
           * (i.logit(X1 %*% par[idx$beta])) * (i.logit(X1 %*% par[idx$theta]))  * (p^Y)*((1-p)^(1-Y)) * ((1-p.prime)^A)*(p.prime^(1-A))

        +  (1/(par[idx$sigma]*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[idx$gamma] + 0*par[idx$gamma_z] ))^2/(2*(par[idx$sigma]^2)) )
           * (1-(i.logit(X1 %*% par[idx$beta]))) * (i.logit(X1 %*% par[idx$theta])) * ((1-p)^Y)*(p^(1-Y))*((1-p.prime)^A)*(p.prime^(1-A))

        +  (1/(par[idx$sigma]*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[idx$gamma] + 1*par[idx$gamma_z] ))^2/(2*(par[idx$sigma]^2)) )
           * (i.logit(X1 %*% par[idx$beta])) * (1 - i.logit(X1 %*% par[idx$theta])) * (1/4)

        +  (1/(par[idx$sigma]*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[idx$gamma] + 0*par[idx$gamma_z] ))^2/(2*(par[idx$sigma]^2)) )
           * (1-(i.logit(X1 %*% par[idx$beta]))) * (1 - i.logit(X1 %*% par[idx$theta])) * (1/4)
  ))
}


# MAXIMIZATION
MLE = optim(par=init,                      # initial values for beta and theta
            fn = log.L.pred,               # function to maximize
            method = "BFGS",               # this method lets set lower bounds (Modified Newton method)
            control = list(maxit=800, fnscale = -1),
            hessian = TRUE)                # calculate Hessian matricce because we will need for confidence intervals

H = MLE$hessian                            # Hessian matrix
Var.hat = diag(-solve(H))                  # Variance as the negative inverse of the Hessian matrix (Dropping covariances)
SE = sqrt(Var.hat)                         # Standard errors


# OUTPUT

  Mlist <- list()
  z = MLE$par / SE
  pv = 2*(1- pnorm(abs(z)))
  z = round(z, d=3)
  pv = round(pv, d=3)


  Mlist[[1]] <- formula
  gamma_idx <- c(idx$gamma, idx$gamma_z)
  theta_sigma_idx <- c(idx$theta, idx$sigma)
  Mlist[[2]] <- t(rbind(MLE$par[gamma_idx],
                             SE[gamma_idx],
                              z[gamma_idx],
                             pv[gamma_idx]))
  Mlist[[3]] <- t(rbind(MLE$par[idx$beta], SE[idx$beta],
                               z[idx$beta], pv[idx$beta]))
  Mlist[[4]] <- t(rbind(MLE$par[theta_sigma_idx],
                             SE[theta_sigma_idx],
                              z[theta_sigma_idx],
                             pv[theta_sigma_idx]))
  Mlist[[2]] <- round(Mlist[[2]], d=4)
  Mlist[[3]] <- round(Mlist[[3]], d=4)
  Mlist[[4]] <- round(Mlist[[4]], d=4)

  colnames(Mlist[[2]]) <- c("Estimate", "Std. Error", "Z score", "Pr(>|z|)")
  colnames(Mlist[[3]]) <- c("Estimate", "Std. Error", "Z score", "Pr(>|z|)")
  colnames(Mlist[[4]]) <- c("Estimate", "Std. Error", "Z score", "Pr(>|z|)")
  Mlist[[5]] <- -solve(H) # Estimated Variance-Covariance Matrix

  varnam <- colnames(X1)
  rownames(Mlist[[2]]) <- c(varnam, crosswise_info$name)
  rownames(Mlist[[3]]) <- varnam
  rownames(Mlist[[4]]) <- c(varnam, "sigma")

names(Mlist) <- c("Call", "Coefficients", "AuxiliaryCoef", "AuxiliaryCoef2", "VCV")

return(Mlist)
}

