#' @title bc.est
#'
#' @description  \code{bc.est} is used to apply a bias-corrected crosswise estimator to survey data.
#'
#' @param Y A vector of binary responses in the crosswise question
#' @param A A vector of binary responses in the anchor question
#' @param p An auxiliary probability for the crosswise question (scalar)
#' @param p.prime An auxiliary probability for the anchor question (scalar)
#' @param w An optional vector specifying sample weights in data (default is set to 1)
#' @param data A dataset containing information from the crosswise model (Y, A, w)
#'
#' @return A list containing main results ($Results) and related statistics ($Stats).
#' @examples
#' bc.est(Y=cross, A=anchor, p=0.15, p.prime=0.15, data=cmdata)
#'
#' #> $Results
#' #>                Point Est. Std. Error Est 95%CI(Lower) 95%CI(Upper)
#' #> Naive Crosswise  0.1950000     0.01444624   0.16668537    0.2233146
#' #> Bias-Corrected   0.1053604     0.02080597   0.06486523    0.1394343
#' #>
#' #> $Stats
#' #> Attentive Rate Est. Sample Size
#' #>           0.7728571        2000
#' @export
#' @importFrom tidyverse

bc.est <- function(Y, A, p, p.prime, w=NULL, data){

Yquo <- enquo(Y)        # QUoting variable name for Y
Aquo <- enquo(A)        # Quoting variable name for A
if(is.null(w)) w <- 1
# data <- data %>% mutate(weight=w)
# data <- data %>% dplyr::select(!!Yquo, !!Aquo, weight)
data <- data[complete.cases(data),] # Just in case (but weight needs to be calculated after dropping NAs)

N = dim(data)[1]        # Number of obs
data <- data %>% as_tibble()
Y = data %>% dplyr::select(!!Yquo) %>% pull()
A = data %>% dplyr::select(!!Aquo) %>% pull()
#weight = data %>% dplyr::select(weight) %>% pull()


# NAIVE CROSSWISE MODEL
  lambda.hat = sum(w*Y)/N                      # Weighted proportion of YESYES or NONO

  pi.hat.naive = (lambda.hat+p-1)/(2*p-1)
  pi.hat.naive = min(1, max(pi.hat.naive, 0))  # Logical bound restrain

  pi.hat.naive.var = (pi.hat.naive*(1-pi.hat.naive))/(N-1) +
    (p*(1-p))/((N-1)*((2*p-1)^2))
  pi.hat.naive.sd = sqrt(pi.hat.naive.var)


  naive.low  = pi.hat.naive-1.96*pi.hat.naive.sd; naive.low  = ifelse(naive.low > 0,  naive.low, 0)
  naive.high = pi.hat.naive+1.96*pi.hat.naive.sd; naive.high = ifelse(naive.high < 1, naive.high, 1)


# BIAS CORRECTED CROSSWISE MODEL
  gamma.hat = (sum(w*A)/N-0.5)/(0.5-p2) # Estimated level of inattentiveness
  Bias.hat = (1/2)*((lambda.hat-0.5)/(p-0.5)) - (1/(2*gamma.hat))*((lambda.hat-0.5)/(p-0.5))
  Bias.hat

  pi.hat.bc = pi.hat.naive - Bias.hat    # Bias Correction
  pi.hat.bc = min(1, max(pi.hat.bc, 0))  # Logical bound restrain


# BOOTSTRAPPING (200 TIMES)
 set.seed(123456)
 bs <- NA
  for(i in 1:200){
    index <- sample(1:nrow(data), size=dim(data), replace=TRUE)   # RESAMPLE WITH REPLACEMENT
    bs.dat <- data[index,]                                        # BOOTSTRAPPED DATA
    N.bs = dim(data)

    Y.bs = bs.dat %>% dplyr::select(!!Yquo) %>% pull()  # First column must be Y
    A.bs = bs.dat %>% dplyr::select(!!Aquo) %>% pull()  # Second column must be A
#    w.bs = bs.dat %>% dplyr::select(weight) %>% pull()  # weight variable
#    w.bs = w.bs/sum(w.bs)                               # normalize the weight
    w.bs=1
    bs.lambda.hat = sum(w.bs*Y.bs)/N.bs                  # Observed proportion of YESYES or NONO
    bs.pi.hat.naive = (bs.lambda.hat+p-1)/(2*p-1)
    bs.gamma.hat = (sum(w.bs*A.bs)/N.bs-0.5)/(0.5-p2)      # Estimated level of inattentiveness
    bs.bias.hat = (1/2)*((bs.lambda.hat-0.5)/(p-0.5)) - (1/(2*bs.gamma.hat))*((bs.lambda.hat-0.5)/(p-0.5))

    bs[i] = bs.pi.hat.naive - bs.bias.hat         # Bias Correction within Bootstrapping
    bs[i] = min(1, max(bs[i], 0))                 # Logical bound restrain
}


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
