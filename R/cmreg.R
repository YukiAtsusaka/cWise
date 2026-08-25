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
#'
#' @return A list containing the estimated results and related statistics.
#' @examples
#' m <- cmreg(Y ~ female + age, anchor = A, p = 0.1, p.prime = 0.15,
#'            data = cmdata2)
#' m
#' @export


cmreg <- function(formula, anchor, p, p.prime, data){

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

  init <- rep(0.01, idx$npar) # Initial values for optim

  # LOG-LIKELIHOOD FUNCTION
  log.L <- function(par) {
    sum(   Y*log(            ((2*p-1)*(i.logit(X1 %*% par[idx$beta])) + (0.5-p) )*i.logit(X1 %*% par[idx$theta]) + 0.5)
           + (1-Y)*log( 1 - (((2*p-1)*(i.logit(X1 %*% par[idx$beta])) + (0.5-p) )*i.logit(X1 %*% par[idx$theta]) + 0.5))
           + A*log(        (0.5-p.prime)*i.logit(X1 %*% par[idx$theta]) + 0.5 )
           + (1-A)*log(1- ((0.5-p.prime)*i.logit(X1 %*% par[idx$theta]) + 0.5 ))
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

  z = MLE$par / SE
  pv = 2*(1- pnorm(abs(z)))
  z = round(z, d=3)
  pv = round(pv, d=3)


  Mlist[[1]] <- formula
  Mlist[[2]] <- t(rbind(MLE$par[idx$beta], SE[idx$beta], z[idx$beta], pv[idx$beta]))
  Mlist[[3]] <- t(rbind(MLE$par[idx$theta], SE[idx$theta],
                        z[idx$theta], pv[idx$theta]))
  Mlist[[2]] <- round(Mlist[[2]], d=4)
  Mlist[[3]] <- round(Mlist[[3]], d=4)
  Mlist[[4]] <- -solve(H) # Estimated Variance-Covariance Matrix
  colnames(Mlist[[2]]) <- c("Estimate", "Std. Error", "z score", "Pr(>|z|)")
  colnames(Mlist[[3]]) <- c("Estimate", "Std. Error", "z score", "Pr(>|z|)")

  varnam <- colnames(X1)
  rownames(Mlist[[2]]) <- varnam
  rownames(Mlist[[3]]) <- varnam

  names(Mlist) <- c("Call", "Coefficients", "AuxiliaryCoef", "VCV")


  return(Mlist)

}

