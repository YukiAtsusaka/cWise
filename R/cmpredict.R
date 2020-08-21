#' @title cmpredict
#'
#' @description Perform a post-estimation prediction with uncertainty quantification via parametric bootstrap
#'
#' @param out An output of \code{cmreg}
#' @param zval A value for the main explanatory variable of interest (the first listed variable in cmreg_out$Call)
#' @param typical A vector of typical values for other covariates
#'
#' @return A vector of predicted probabilities given the input covariatevalues
#' @examples
#' m <- cmpredict(m, zval=0, typical=30)
#' m
#' @export


cmpredict <- function(out, zval, typical){

# GRAB COEFFICIENTS
  coef.beta  = out$Coefficients[,1]
  coef.theta = out$AuxiliaryCoef[,1]
  coefs = c(coef.beta, coef.theta)
  vcovs = out$VCV

# TYPICAL VALUE MATRIX
  typ.vec = cbind(1, zval, typical)

# PARAMETRIC BOOTSTRAP
  k = 1 + length(typical) + 1
  set.seed(20200730)
  coef.sim <- mvtnorm::rmvnorm(n=10000, mean=coefs, sigma=vcovs) # from mvtnorm
  coef.sim <- coef.sim[,1:k]

  lin.agg <- as.matrix(typ.vec) %*% t(coef.sim) # Linear aggregator
  pi.sim = i.logit(lin.agg) # Inverse logit

  return(pi.sim)
}
