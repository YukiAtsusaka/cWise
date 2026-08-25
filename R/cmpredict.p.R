#' @title cmpredict.p
#'
#' @description Perform a post-estimation prediction with uncertainty quantification via parametric bootstrap
#'
#' @param out An output of \code{cmreg.p}
#' @param zval A value for the explanatory variable of interest (the first listed variable in cmreg_out$Call)
#' @param typical A vector of typical values for other covariates
#'
#' @return A vector of predicted probabilities given the input covariatevalues
#' @examples
#' pr2 <- cmpredict2(out=m2, typical=c(1,30))
#' pr2
#' @export
#' @importFrom mvtnorm "rmvnorm"

cmpredict.p <- function(out, typical){

# GRAB COEFFICIENTS
  idx <- cm_par_index(nrow(out$AuxiliaryCoef), model = "predictor")
  gamma_idx = c(idx$gamma, idx$gamma_z)
  coef.gamma = out$Coefficients[,1]
  vcovs = out$VCV[gamma_idx, gamma_idx, drop = FALSE]


# PARAMETRIC BOOTSTRAP
  set.seed(20200730)
  coef.sim <- rmvnorm(n=10000, mean=coef.gamma, sigma=vcovs) # from mvtnorm


# TYPICAL VALUE MATRIX
  typ.vec = rbind(c(1, typical, 0),
                  c(1, typical, 1))

  lin.agg <- as.matrix(typ.vec) %*% t(coef.sim) # Linear aggregator

  return(lin.agg)
}

