#' @title cmpredict
#'
#' @description Perform a post-estimation prediction with uncertainty quantification via parametric bootstrap
#'
#' @param formula explain
#' @param p explain
#' @param p.prime explain
#' @param data explain
#' @param init explain
#'
#' @return ggplot object
#' @examples
#' m <-  cmreg(Y~female+age+A, p=0.1, p.prime=0.15, data=cmdata2)
#' m
#' @export
#' @importFrom dplyr


cmpredict <- function(cmreg_obj, typical, zval){

# GRAB COEFFICIENTS
  coef.beta = cmreg_obj$Coefficients[,1]
  coef.theta = cmreg_obj$AuxiliaryCoef[,1]
  coefs = c(coef.beta, coef.theta)
  vcovs = cmreg_obj$VCV

# TYPICAL VALUE MATRIX
  typ.vec = cbind(1, zval, typical)


# PARAMETRIC BOOTSTRAP
  k = 1 + length(typical) + 1
  set.seed(20200730)
  coef.sim <- rmvnorm(n=10000, mean=coefs, sigma=vcovs) # from mvtnorm
  coef.sim <- coef.sim[,1:k]

  lin.agg <- as.matrix(typ.vec) %*% t(coef.sim) # Linear aggregator
  pi.sim = i.logit(lin.agg) # Inverse logit

  return(pi.sim)
}

