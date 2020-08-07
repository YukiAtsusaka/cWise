#' @title cmreg
#'
#' @description \code{cmreg} is used to run a regression with the latent sensitive trait as an outcome.
#'
#' @param formula an object of class "formula":
#' a symbolic description of the model to be fitted.
#' Ex. Crosswise response ~ Covariates + Anchor response.
#' @param p an auxiliary probability for the crosswise question.
#' @param p.prime an auxiliary probability for the anchor question.
#' @param data a data frame containing information from the crosswise model and covariates.
#'
#' @return A list containing the estimated results and related statistics.
#' @examples
#' m <-  cmreg(Y~female+age+A, p=0.1, p.prime=0.15, data=cmdata2)
#' m
#' @export


cmreg <- function(formula, p, p.prime, data){

  i.logit <- function(XB){ exp(XB)/(1 + exp(XB))}

  df <- model.frame(formula, data, na.action = na.omit)
  X0 <- model.matrix.default(formula, df) # Data matrix including A
  X1 <- X0[, -dim(X0)[2]]                 # Matrix with 1 and predictors
  A <- df[, dim(X0)[2]]
  Y <- df[,1]
  k <- dim(X1)[2]      # Number of beta parameters
  k.t <- k + 1         # Start of theta parameters
  k.t.e <- 2*k         # End of theta parameters

  init <- rep(0.01, k.t.e) # Initial values for optim

  # LOG-LIKELIHOOD FUNCTION
  log.L <- function(par) {
    sum(   Y*log(            ((2*p-1)*(i.logit(X1 %*% par[1:k])) + (0.5-p) )*i.logit(X1 %*% par[k.t:k.t.e]) + 0.5)
           + (1-Y)*log( 1 - (((2*p-1)*(i.logit(X1 %*% par[1:k])) + (0.5-p) )*i.logit(X1 %*% par[k.t:k.t.e]) + 0.5))
           + A*log(        (0.5-p.prime)*i.logit(X1 %*% par[k.t:k.t.e]) + 0.5 )
           + (1-A)*log(1- ((0.5-p.prime)*i.logit(X1 %*% par[k.t:k.t.e]) + 0.5 ))
    )
  }


  # MAXIMIZATION
  MLE = optim(par=init,                      # initial values for beta and theta
              fn = log.L,                    # function to maximize
              method = "BFGS",               # this method lets set lower bounds (Modified Newton method)
              control = list(maxit=800, fnscale = -1),  # maximize the function
              hessian = TRUE)                # calculate Hessian matricce because we will need for confidence intervals

  H = MLE$hessian                            # Hessian matrix
  Var.hat = diag(-solve(H))                  # Variance as the negative inverse of the Hessian matrix (Dropping covariances)
  SE = sqrt(Var.hat)                         # Standard errors


# OUTPUT

  Mlist <- list()

  n.var = dim(df)[2] - 1
  z = MLE$par / SE
  pv = 2*(1- pnorm(abs(z)))
  z = round(z, d=3)
  pv = round(pv, d=3)


  Mlist[[1]] <- formula
  Mlist[[2]] <- t(rbind(MLE$par[1:n.var], SE[1:n.var], z[1:n.var], pv[1:n.var]))
  Mlist[[3]] <- t(rbind(MLE$par[(n.var+1):(2*n.var)], SE[(n.var+1):(2*n.var)],
                        z[(n.var+1):(2*n.var)], pv[(n.var+1):(2*n.var)]))
  Mlist[[2]] <- round(Mlist[[2]], d=4)
  Mlist[[3]] <- round(Mlist[[3]], d=4)
  Mlist[[4]] <- -solve(H) # Estimated Variance-Covariance Matrix
  colnames(Mlist[[2]]) <- c("Estimate", "Std. Error", "z score", "Pr(>|z|)")
  colnames(Mlist[[3]]) <- c("Estimate", "Std. Error", "z score", "Pr(>|z|)")

  varnam <- c("(intercept)", colnames(df)[2:n.var])
  rownames(Mlist[[2]]) <- varnam
  rownames(Mlist[[3]]) <- varnam

  names(Mlist) <- c("Call", "Coefficients", "AuxiliaryCoef", "VCV")


  return(Mlist)

}

