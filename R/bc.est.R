#' @title bc.est
#'
#' @description Apply a bias-corrected crosswise estimator to survey data
#'
#' @param Y explain
#' @param A explain
#' @param p explain
#' @param p.prime explain
#'
#' @return ggplot object
#' @examples
#' sensitivity <- cmBound(p=0.25, lambda.hat=0.6385, N=310, dq=0.073)
#' @export
#' @importFrom dplyr


bc.est <- function(Y, A, p, p2, data){

Yquo <- enquo(Y)        # QUoting variable name for Y
Aquo <- enquo(A)        # Quoting variable name for A
data <- data %>% dplyr::select(!!Yquo, !!Aquo)
data <- data[complete.cases(data),]

N = dim(data)[1]        # Number of obs
data <- data %>% as_tibble()
Y = data %>% dplyr::select(!!Yquo) %>% pull()
A = data %>% dplyr::select(!!Aquo) %>% pull()


# NAIVE CROSSWISE MODEL
  lambda.hat = mean(Y)                         # Observed proportion of YESYES or NONO
  pi.hat.naive = (lambda.hat+p-1)/(2*p-1)
  pi.hat.naive = min(1, max(pi.hat.naive, 0))  # Logical bound restrain

  pi.hat.naive.var = (pi.hat.naive*(1-pi.hat.naive))/(N-1) +
    (p*(1-p))/((N-1)*((2*p-1)^2))
  pi.hat.naive.sd = sqrt(pi.hat.naive.var)


  naive.low  = pi.hat.naive-1.96*pi.hat.naive.sd; naive.low  = ifelse(naive.low > 0,  naive.low, 0)
  naive.high = pi.hat.naive+1.96*pi.hat.naive.sd; naive.high = ifelse(naive.high < 1, naive.high, 1)



# BIAS CORRECTED CROSSWISE MODEL
  gamma.hat = (mean(A)-0.5)/(0.5-p2) # Estimated level of inattentiveness
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

    Y.bs = bs.dat %>% dplyr::select(!!Yquo) %>% pull()  # First column must be Y
    A.bs = bs.dat %>% dplyr::select(!!Aquo) %>% pull()  # Second column must be A

    bs.lambda.hat = mean(Y.bs)                    # Observed proportion of YESYES or NONO
    bs.pi.hat.naive = (bs.lambda.hat+p-1)/(2*p-1)
    bs.gamma.hat = (mean(A.bs)-0.5)/(0.5-p2) # Estimated level of inattentiveness
    bs.bias.hat = (1/2)*((bs.lambda.hat-0.5)/(p-0.5)) - (1/(2*bs.gamma.hat))*((bs.lambda.hat-0.5)/(p-0.5))

    bs[i] = bs.pi.hat.naive - bs.bias.hat # Bias Correction within Bootstrapping
    bs[i] = min(1, max(bs[i], 0))  # Logical bound restrain

  }

  pi.hat.bc.var = var(bs)
  pi.hat.bc.sd = sqrt(pi.hat.bc.var)

  bc.low  = quantile(bs, prob=0.025); bc.low  = ifelse(bc.low > 0,  bc.low, 0)
  bc.high = quantile(bs, prob=0.975); bc.high = ifelse(bc.high < 1, bc.high, 1)


result <-   rbind(cbind(pi.hat.naive, pi.hat.naive.sd, naive.low, naive.high, N),
                  cbind(pi.hat.bc, pi.hat.bc.sd, bc.low, bc.high, N))

return(result)
}
