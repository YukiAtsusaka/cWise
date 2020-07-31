#' @title cmpredict2
#'
#' @description Perform a post-estimation prediction with uncertainty quantification via parametric bootstrap
#'
#' @param cmreg_out An output of \code{cmreg} or \code{cmreg.p}
#' @param zval A value for the explanatory variable of interest (the first listed variable in cmreg_out$Call)
#' @param typical A vector of typical values for other covariates
#'
#' @return A vector of predicted probabilities given the input covariatevalues
#' @examples
#' pr2 <- cmpredict2(out=m2, typical=c(1,30))
#' pr2
#' @export
#' @importFrom dplyr


cmpredict2 <- function(out, typical){

# GRAB COEFFICIENTS
  coef.gamma = out$Coefficients[,1]
  coef.beta  = out$AuxiliaryCoef[,1]
  coef.theta = out$AuxiliaryCoef2[,1]
  coefs = c(coef.gamma, coef.beta, coef.theta)
  vcovs = out$VCV


# PARAMETRIC BOOTSTRAP
  k = 1 + length(typical) + 1
  set.seed(20200730)
  coef.sim <- rmvnorm(n=10000, mean=coefs, sigma=vcovs) # from mvtnorm
  coef.sim <- coef.sim[,1:k]


  # TYPICAL VALUE MATRIX
  typ.vec = rbind(c(1, typical, 0),
                  c(1, typical, 1))


  lin.agg <- as.matrix(typ.vec) %*% t(coef.sim) # Linear aggregator

  return(lin.agg)
}

#
# pr2 <- cmpredict2(m2, typical=c(1,30))
# #
# par(mfrow=c(1,2))
# hist(pr2[1,], main="No Sensitive Trait", xlab="Outcome Value", breaks=40)
# hist(pr2[2,], main="With Sensitive Trait", xlab="Outcome Value", breaks=40)


