#' @title cmreg.p
#'
#' @description Run a regression with the latent sensitive trait as a predictor
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


cmreg.p <- function(formula, p, p2, data){

  i.logit <- function(XB){ exp(XB)/(1 + exp(XB))}
  df <- model.frame(formula, data, na.action = na.omit)
  X0 <- model.matrix.default(formula, df) # Data matrix including A
  X1 <- X0[, 1:(dim(X0)[2]-2)]            # Matrix with 1 and predictors

  A <- df[, dim(X0)[2]]
  Y <- df[,dim(X0)[2]-1]
  V <- df[,1]          # Outcome variable
  p <- p               # Randomization probability in Crosswise Q
  p2 <- p2             # Randomization probability in Anchor Q

  k <- dim(X1)[2]      # Number of beta parameters ( #covariate + 1)
  k.t <- k + 1         # Start of theta parameters
  k.t.e <- 2*k         # End of theta parameters
  k.g <- k.t.e + 1     # Start of gamma parameters ( #coveriate + 2)
  k.g.e <- k.g + k + 1 # End of gamma parameters
  k.g.e_part <- k.g.e - 2

# par[k.g.e] = sigma parameter
# par[k.g.e-1]  = gamma.cm (coef on Z)

  init <- c(rep(0.01,k.g.e-1),1)

#---------------------------------------------------------------------------------------#
# PARAMTERS TO ESTIMATE (11 parameters, if two covariates)
# beta0, beta1, beta2,
#  theta0, theta1, theta2,
#  gamma0, gamma1, gamma2, gamma.cm (Z), sigma
#---------------------------------------------------------------------------------------#

# LOG-LIKELIHOOD FUNCTION
log.L.pred <- function(par) {
  sum(log( (1/(par[k.g.e]*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[k.g:k.g.e_part] + 1*par[k.g.e-1] ))^2/(2*(par[k.g.e]^2)) )
           * (i.logit(X1 %*% par[1:k])) * (i.logit(X1 %*% par[k.t:k.t.e]))  * (p^Y)*((1-p)^(1-Y)) * ((1-p2)^A)*(p2^(1-A))

        +  (1/(par[k.g.e]*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[k.g:k.g.e_part] + 0*par[k.g.e-1] ))^2/(2*(par[k.g.e]^2)) )
           * (1-(i.logit(X1 %*% par[1:k]))) * (i.logit(X1 %*% par[k.t:k.t.e])) * ((1-p)^Y)*(p^(1-Y))*((1-p2)^A)*(p2^(1-A))

        +  (1/(par[k.g.e]*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[k.g:k.g.e_part] + 1*par[k.g.e-1] ))^2/(2*(par[k.g.e]^2)) )
           * (i.logit(X1 %*% par[1:k])) * (1 - i.logit(X1 %*% par[k.t:k.t.e])) * (1/4)

        +  (1/(par[k.g.e]*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[k.g:k.g.e_part] + 0*par[k.g.e-1] ))^2/(2*(par[k.g.e]^2)) )
           * (1-(i.logit(X1 %*% par[1:k]))) * (1 - i.logit(X1 %*% par[k.t:k.t.e])) * (1/4)
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


Mlist <- list()
Mlist[[1]] <- MLE$par
Mlist[[2]] <- SE

return(Mlist)
# Next: make pretty summary
}

