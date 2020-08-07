#' @title bc.est
#'
#' @description  \code{bc.est} is used to apply a bias-corrected crosswise estimator to survey data.
#'
#' @param Y a vector of binary responses in the crosswise question (Y=1 if TRUE-TRUE or FALSE-FALSE, Y=0 otherwise).
#' @param A a vector of binary responses in the anchor question (A=1 if TRUE-TRUE or FALSE-FALSE, A=0 otherwise).
#' @param p an auxiliary probability for the crosswise question.
#' @param p.prime an auxiliary probability for the anchor question.
#' @param weight an optional vector specifying sample weights in data.
#' @param data a data frame containing information from the crosswise model (Y, A, weight).
#'
#' @return A list containing main results ($Results) and related statistics ($Stats).
#' @examples
#' bc.est(Y=Y, A=A, p=0.15, p.prime=0.15, data=cmdata)
#' bc.est(Y=Y, A=A, weight=weight, p=0.15, p.prime=0.15, data=cmdata)
#' @export
#' @importFrom dplyr "%>%" "select" "pull"
#'             tibble "as_tibble"

bc.est <- function(Y, A, p, p.prime, weight, data){

Yquo <- dplyr::enquo(Y)        # QUoting variable name for Y
Aquo <- dplyr::enquo(A)        # Quoting variable name for A
data <- data[complete.cases(data),] # Just in case (but weight needs to be calculated after dropping NAs)

N = dim(data)[1]        # Number of obs
data <- data %>% as_tibble()
Y = data %>% dplyr::select(!!Yquo) %>% pull()
A = data %>% dplyr::select(!!Aquo) %>% pull()

if(missing(weight)){
  weight <- rep(1, N)    # If weight is not specified
}else{
Wquo <- dplyr::enquo(weight)   # Quoting variable name for weight
weight = data %>% dplyr::select(!!Wquo) %>% pull()
}

# NAIVE CROSSWISE MODEL
  lambda.hat = sum(weight*Y)/ sum(weight)                      # Weighted proportion of YESYES or NONO

  pi.hat.naive = (lambda.hat+p-1)/(2*p-1)
  pi.hat.naive = min(1, max(pi.hat.naive, 0))  # Logical bound restrain

  pi.hat.naive.var = (pi.hat.naive*(1-pi.hat.naive))/(N-1) +
    (p*(1-p))/((N-1)*((2*p-1)^2))
  pi.hat.naive.sd = sqrt(pi.hat.naive.var)


  naive.low  = pi.hat.naive-1.96*pi.hat.naive.sd; naive.low  = ifelse(naive.low > 0,  naive.low, 0)
  naive.high = pi.hat.naive+1.96*pi.hat.naive.sd; naive.high = ifelse(naive.high < 1, naive.high, 1)


# BIAS CORRECTED CROSSWISE MODEL
  gamma.hat = (sum(weight*A) / sum(weight)-0.5)/(0.5-p2) # Estimated level of inattentiveness
  Bias.hat = (1/2)*((lambda.hat-0.5)/(p-0.5)) - (1/(2*gamma.hat))*((lambda.hat-0.5)/(p-0.5))
  Bias.hat

  pi.hat.bc = pi.hat.naive - Bias.hat    # Bias Correction
  pi.hat.bc = min(1, max(pi.hat.bc, 0))  # Logical bound restrain


# BOOTSTRAPPING FOR UNCERTAINTY ESTIMATE (200 TIMES)
 set.seed(123456)
 bs <- NA
  for(i in 1:200){
    index <- sample(1:nrow(data), size=dim(data)[1], replace=TRUE)   # RESAMPLE WITH REPLACEMENT
    bs.dat <- data[index,]                                        # BOOTSTRAPPED DATA
    N.bs = dim(data)[1]

    Y.bs = bs.dat %>% dplyr::select(!!Yquo) %>% pull()  # First column must be Y
    A.bs = bs.dat %>% dplyr::select(!!Aquo) %>% pull()  # Second column must be A

    w.bs = weight[index] # Carefully see if this is working
# if(missing(weight)){
# w.bs <- rep(1, N.bs)    # If weight is not specified
# }else{
# Wquo <- enquo(weight)   # Quoting variable name for weight
# w.bs = bs.dat %>% dplyr::select(!!Wquo) %>% pull()
# }

    bs.lambda.hat = sum(w.bs*Y.bs)/sum(w.bs)                  # Observed proportion of YESYES or NONO
    bs.pi.hat.naive = (bs.lambda.hat+p-1)/(2*p-1)
    bs.gamma.hat = (sum(w.bs*A.bs)/sum(w.bs)-0.5)/(0.5-p2)      # Estimated level of inattentiveness
    bs.bias.hat = (1/2)*((bs.lambda.hat-0.5)/(p-0.5)) - (1/(2*bs.gamma.hat))*((bs.lambda.hat-0.5)/(p-0.5))

    bs[i] = bs.pi.hat.naive - bs.bias.hat         # Bias Correction within Bootstrapping
    bs[i] = min(1, max(bs[i], 0))                 # Logical bound restrain
} # END OF BOOTSTRAPPING


  pi.hat.bc.var = var(bs)
  pi.hat.bc.sd = sqrt(pi.hat.bc.var)

  bc.low  = quantile(bs, prob=0.025); bc.low  = ifelse(bc.low > 0,  bc.low, 0)
  bc.high = quantile(bs, prob=0.975); bc.high = ifelse(bc.high < 1, bc.high, 1)


# Output
  result <- list()

  result[[1]] <-   rbind(cbind(pi.hat.naive, pi.hat.naive.sd, naive.low, naive.high),
                         cbind(pi.hat.bc, pi.hat.bc.sd, bc.low, bc.high))
  result[[1]] <- round(result[[1]], d=4)
  colnames(result[[1]]) <- c("Estimate", "Std. Error", "95%CI(Low)", "95%CI(Up)")
  rownames(result[[1]]) <- c("Naive Crosswise", "Bias-Corrected")

  result[[2]] <- cbind(gamma.hat, N)
  result[[2]] <- round(result[[2]], d=4)
  colnames(result[[2]]) <- c("Attentive Rate", "Sample Size")
  rownames(result[[2]]) <- ""
  names(result) <- c("Results", "Stats")

return(result)
}
