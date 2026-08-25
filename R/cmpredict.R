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
#' @importFrom mvtnorm "rmvnorm"


cmpredict <- function(out, zval, typical){

i.logit <- function(XB){ exp(XB)/(1 + exp(XB))}

# GRAB COEFFICIENTS
  idx <- cm_par_index(nrow(out$Coefficients), model = "outcome")
  coef.beta  = out$Coefficients[,1]
  vcovs = out$VCV[idx$beta, idx$beta, drop = FALSE]

# TYPICAL VALUE MATRIX
  typ.vec = cbind(1, zval, typical)

# PARAMETRIC BOOTSTRAP
  set.seed(20200730)
  coef.sim <- mvtnorm::rmvnorm(n=10000, mean=coef.beta, sigma=vcovs) # from mvtnorm

  lin.agg <- as.matrix(typ.vec) %*% t(coef.sim) # Linear aggregator
  pi.sim = i.logit(lin.agg) # Inverse logit

  return(pi.sim)
}
