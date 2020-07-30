#' @title cmreg.p
#'
#' @description Run a regression with the latent sensitive trait as a predictor
#'
#' @param formula explain
#' @param p explain
#' @param p.prime explain
#' @param data explain
#'
#' @return ggplot object
#' @examples
#' m2 <- cmreg.p(V~age+female+Y+A, p=0.1, p.prime=0.15, data=cmdata3)
#' m2
#' @export
#' @importFrom dplyr


cmreg.p <- function(formula, p, p.prime, data){
  options(warn=-1)

  i.logit <- function(XB){ exp(XB)/(1 + exp(XB))}
  df <- model.frame(formula, data, na.action = na.omit)
  X0 <- model.matrix.default(formula, df) # Data matrix including A
  X1 <- X0[, 1:(dim(X0)[2]-2)]            # Matrix with 1 and predictors

  A <- df[, dim(X0)[2]]
  Y <- df[,dim(X0)[2]-1]
  V <- df[,1]          # Outcome variable

  k <- dim(X1)[2]      # Number of beta parameters ( #covariate + 1)
  k.t <- k + 1         # Start of theta parameters
  k.t.e <- 2*k         # End of theta parameters
  k.g <- k.t.e + 1     # Start of gamma parameters ( #coveriate + 2)
  k.g.e <- k.g + k + 1 # End of gamma parameters
  k.g.e_part <- k.g.e - 2

  init <- c(rep(0.01,k.g.e-1),1) # Initial values for optimization

#---------------------------------------------------------------------------------------#
# PARAMTERS TO ESTIMATE (11 parameters, if two covariates)
# beta0, beta1, beta2,
# theta0, theta1, theta2,
# gamma0, gamma1, gamma2, gamma.cm (Z), sigma
#---------------------------------------------------------------------------------------#

# LOG-LIKELIHOOD FUNCTION
log.L.pred <- function(par) {
  sum(log( (1/(par[k.g.e]*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[k.g:k.g.e_part] + 1*par[k.g.e-1] ))^2/(2*(par[k.g.e]^2)) )
           * (i.logit(X1 %*% par[1:k])) * (i.logit(X1 %*% par[k.t:k.t.e]))  * (p^Y)*((1-p)^(1-Y)) * ((1-p.prime)^A)*(p.prime^(1-A))

        +  (1/(par[k.g.e]*sqrt(2*pi)) )*exp( -(V - (X1 %*% par[k.g:k.g.e_part] + 0*par[k.g.e-1] ))^2/(2*(par[k.g.e]^2)) )
           * (1-(i.logit(X1 %*% par[1:k]))) * (i.logit(X1 %*% par[k.t:k.t.e])) * ((1-p)^Y)*(p^(1-Y))*((1-p.prime)^A)*(p.prime^(1-A))

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


# OUTPUT

  Mlist <- list()
  n.var = dim(df)[2] - 2

  Mlist[[1]] <- formula
  Mlist[[2]] <- t(rbind(MLE$par[(2*n.var+1):(3*n.var+1)], SE[(2*n.var+1):(3*n.var+1)]))
  Mlist[[3]] <- t(rbind(MLE$par[1:n.var], SE[1:n.var]))
  Mlist[[4]] <- t(rbind(MLE$par[(n.var+1):(2*n.var)], SE[(n.var+1):(2*n.var)]))
  Mlist[[2]] <- round(Mlist[[2]], d=4)
  Mlist[[3]] <- round(Mlist[[3]], d=4)
  Mlist[[4]] <- round(Mlist[[4]], d=4)

  colnames(Mlist[[2]]) <- c("Estimate", "Std. Error")
  colnames(Mlist[[3]]) <- c("Estimate", "Std. Error")
  colnames(Mlist[[4]]) <- c("Estimate", "Std. Error")

  varnam   <- c("(intercept)", colnames(df)[2:(n.var+1)])
  varnam.s <- c("(intercept)", colnames(df)[2:(n.var)])
  rownames(Mlist[[2]]) <- varnam
  rownames(Mlist[[3]]) <- varnam.s
  rownames(Mlist[[4]]) <- varnam.s

names(Mlist) <- c("Call", "Coefficients", "AuxiliaryCoef", "AuxiliaryCoef2")

return(Mlist)
}

