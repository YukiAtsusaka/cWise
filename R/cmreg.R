#' @title cmreg
#'
#' @description Run a regression with the latent sensitive trait as an outcome
#'
#' @param formula explain
#' @param p explain
#' @param p.prime explain
#' @param data explain
#' @param init explain
#'
#' @return ggplot object
#' @examples
#' sensitivity <- cmBound(p=0.25, lambda.hat=0.6385, N=310, dq=0.073)
#' @export
#' @importFrom dplyr


cmreg <- function(formula, p, p2, data, init){
i.logit <- function(XB){ exp(XB)/(1 + exp(XB))}

  df <- model.frame(formula, data, na.action = na.omit)
  X0 <- model.matrix.default(formula, df) # Data matrix including A
  X1 <- X0[, -dim(X0)[2]]                 # Matrix with 1 and predictors
  A <- df[, dim(X0)[2]]
  Y <- df[,1]
  p <- p
  p2 <- p2
  k <- dim(X1)[2]      # Number of beta parameters
  k.t <- k + 1         # Start of theta parameters
  k.t.e <- 2*k         # End of theta parameters

  init <- init         # Initial values for optim

  # LOG-LIKELIHOOD FUNCTION
  log.L <- function(par) {
    sum(   Y*log(            ((2*p-1)*(i.logit(X1 %*% par[1:k])) + (0.5-p) )*i.logit(X1 %*% par[k.t:k.t.e]) + 0.5)
           + (1-Y)*log( 1 - (((2*p-1)*(i.logit(X1 %*% par[1:k])) + (0.5-p) )*i.logit(X1 %*% par[k.t:k.t.e]) + 0.5))
           + A*log(        (0.5-p2)*i.logit(X1 %*% par[k.t:k.t.e]) + 0.5 )
           + (1-A)*log(1- ((0.5-p2)*i.logit(X1 %*% par[k.t:k.t.e]) + 0.5 ))
    )
  }


  # MAXIMIZATION
  MLE = optim(par=init,                   # initial values for beta and theta
              fn = log.L,                    # function to maximize
              method = "BFGS",               # this method lets set lower bounds (Modified Newton method)
              control = list(maxit=800, fnscale = -1),  # maximize the function
              hessian = TRUE)                # calculate Hessian matricce because we will need for confidence intervals

  H = MLE$hessian                            # Hessian matrix
  Var.hat = diag(-solve(H))                  # Variance as the negative inverse of the Hessian matrix (Dropping covariances)
  SE = sqrt(Var.hat)                         # Standard errors


  Mlist <- list()
  Mlist[[1]] <- MLE$par
  Mlist[[2]] <- SE

  return(Mlist)

  # Next: make pretty summary
}

