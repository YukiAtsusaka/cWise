#' Simulate Survey Data and Compute Bias-Corrected Estimates
#'
#' Core simulation function that generates crosswise model data with inattentive
#' respondents and computes bias-corrected prevalence estimates with bootstrap
#' confidence intervals. This implements the methodology described in Appendix C5
#' for power analysis of the crosswise model.
#'
#' @param N.sim Integer. Number of Monte Carlo simulations to run. Default is 500.
#' @param sample Integer. Sample size per simulation (number of respondents).
#' @param pi Numeric. True prevalence rate of the sensitive attribute (between 0 and 1).
#' @param p Numeric. Probability for the randomization item in sensitive question (between 0 and 1).
#' @param p.prime Numeric. Probability for the anchor question (non-sensitive, between 0 and 1).
#' @param gamma Numeric. Proportion of attentive respondents (between 0 and 1).
#' @param direct Numeric. Direct questioning estimate for comparison purposes (between 0 and 1).
#'
#' @return A list containing:
#' \describe{
#'   \item{Results}{Named vector with summary statistics including average bias,
#'                  RMSE, and coverage rates for both naive and bias-corrected estimators}
#'   \item{BiasCorrectEst}{Numeric vector of bias-corrected point estimates (sorted)}
#'   \item{BiasCorrectLow}{Numeric vector of lower bounds of 95\% bootstrap CIs (sorted)}
#'   \item{BiasCorrectHigh}{Numeric vector of upper bounds of 95\% bootstrap CIs (sorted)}
#'   \item{EstimatedBias}{Numeric vector of estimated bias values for each simulation}
#'   \item{RelativeLengthCI}{Numeric vector of relative CI lengths (bias-corrected vs naive)}
#' }
#'
#' @details
#' The function implements the crosswise model with bias correction for inattentive
#' respondents. For each simulation:
#' \itemize{
#'   \item Generates attentive/inattentive status based on gamma
#'   \item Simulates responses to sensitive question (crosswise format)
#'   \item Simulates responses to anchor question (for estimating attention rate)
#'   \item Computes naive crosswise estimate
#'   \item Applies bias correction based on estimated inattention
#'   \item Uses bootstrap (500 iterations) to compute 95\% confidence intervals
#' }
#'
#' The bias correction formula is:
#' \deqn{\hat{\pi}_{BC} = \hat{\pi}_{naive} - \hat{Bias}}
#' where the bias is estimated using the anchor question responses.
#'
#' @examples
#' # Basic usage
#' result <- sim.power(
#'   N.sim = 100,
#'   sample = 500,
#'   pi = 0.1,
#'   p = 0.1,
#'   p.prime = 0.1,
#'   gamma = 0.8,
#'   direct = 0.05
#' )
#' print(result$Results)
#'
#' @references
#' Atsusaka and Stevenson (2021). Appendix C5: Sample Size Determination
#' and Parameter Selection.
#'
#' @export
sim.power <- function(N.sim = 500, sample, pi, p, p.prime, gamma, direct) {

  # Helper function for rounding output
  roundy <- function(x) { round(x, digits = 3) }

  N <- sample
  p2 <- p.prime
  pb <- txtProgressBar(min = 1, max = N.sim, initial = 0, style = 3)

  # Storage for simulation results
  naive.cover <- numeric(N.sim)
  bc.cover <- numeric(N.sim)
  naive.pred <- numeric(N.sim)
  bc.pred <- numeric(N.sim)
  bc.high.save <- numeric(N.sim)
  bc.low.save <- numeric(N.sim)
  REbc <- numeric(N.sim)
  est.bias <- numeric(N.sim)
  bc.cover0 <- numeric(N.sim)
  bc.cover.direct <- numeric(N.sim)

  #############################################################
  # Repeat simulation N.sim times with the same parameters
  #############################################################
  for (i in 1:N.sim) {

    setTxtProgressBar(pb, i)

    # Generate attentive/inattentive status
    Attentive <- rbinom(n = N, size = 1, prob = gamma)

    # SENSITIVE QUESTION OF INTEREST (Crosswise format)
    StatementA <- rbinom(n = N, size = 1, prob = pi)           # Sensitive item
    StatementB <- rbinom(n = N, size = 1, prob = p)            # Randomization item
    True.res <- ifelse(StatementA != StatementB, 0, 1)         # True answer (same=1, different=0)

    Obs.res <- rep(NA, N)
    Obs.res[Attentive == 1] <- True.res[Attentive == 1]        # Attentive Rs give true answer
    Obs.res[Attentive == 0] <- rbinom(n = sum(Attentive == 0),
                                      size = 1, prob = 0.5)     # Inattentive Rs answer randomly

    # NON-SENSITIVE AUXILIARY QUESTION (Anchor question for estimating gamma)
    StatementC <- 0                                             # Known prevalence of 0
    StatementD <- rbinom(n = N, size = 1, prob = p2)            # Randomization item
    True.res.prime <- ifelse(StatementC != StatementD, 0, 1)

    Obs.res.prime <- rep(NA, N)
    Obs.res.prime[Attentive == 1] <- True.res.prime[Attentive == 1]
    Obs.res.prime[Attentive == 0] <- rbinom(n = sum(Attentive == 0),
                                            size = 1, prob = 0.5)

    dat <- data.frame(Obs.res, True.res, Attentive, Obs.res.prime, True.res.prime)


    # (2) NAIVE CROSSWISE MODEL ESTIMATOR
    lambda.hat <- mean(Obs.res)                                 # Observed proportion
    pi.hat.naive <- (lambda.hat + p - 1) / (2*p - 1)

    pi.hat.naive.var <- (pi.hat.naive * (1 - pi.hat.naive)) / (N - 1) +
                        (p * (1 - p)) / ((N - 1) * ((2*p - 1)^2))
    pi.hat.naive.sd <- sqrt(pi.hat.naive.var)

    naive.low  <- max(0, pi.hat.naive - 1.96 * pi.hat.naive.sd)
    naive.high <- min(1, pi.hat.naive + 1.96 * pi.hat.naive.sd)


    # (3) BIAS CORRECTED CROSSWISE MODEL ESTIMATOR
    gamma.hat <- (mean(Obs.res.prime) - 0.5) / (0.5 - p2)       # Estimated attention rate
    Bias.hat <- (1/2) * ((lambda.hat - 0.5) / (p - 0.5)) -
                (1 / (2 * gamma.hat)) * ((lambda.hat - 0.5) / (p - 0.5))

    pi.hat.bc <- pi.hat.naive - Bias.hat                        # Apply bias correction
    pi.hat.bc <- max(0, min(1, pi.hat.bc))                      # Bound between [0,1]

    # BOOTSTRAPPING for confidence intervals (500 iterations)
    bs <- numeric(500)
    for (b in 1:500) {
      index <- sample(1:nrow(dat), size = nrow(dat), replace = TRUE)
      bs.dat <- dat[index, ]

      bs.lambda.hat <- mean(bs.dat$Obs.res)
      bs.pi.hat.naive <- (bs.lambda.hat + p - 1) / (2*p - 1)
      bs.gamma.hat <- (mean(bs.dat$Obs.res.prime) - 0.5) / (0.5 - p2)
      bs.bias.hat <- (1/2) * ((bs.lambda.hat - 0.5) / (p - 0.5)) -
                     (1 / (2 * bs.gamma.hat)) * ((bs.lambda.hat - 0.5) / (p - 0.5))

      bs[b] <- bs.pi.hat.naive - bs.bias.hat
      bs[b] <- max(0, min(1, bs[b]))                            # Bound between [0,1]
    } # END OF BOOTSTRAPPING

    bc.low  <- max(0, quantile(bs, prob = 0.025, na.rm = TRUE))
    bc.high <- min(1, quantile(bs, prob = 0.975, na.rm = TRUE))


    # STORE RESULTS from this simulation
    est.bias[i] <- Bias.hat
    naive.pred[i] <- pi.hat.naive
    bc.pred[i] <- pi.hat.bc
    bc.high.save[i] <- bc.high
    bc.low.save[i] <- bc.low
    naive.cover[i] <- (naive.low <= pi & pi <= naive.high)
    bc.cover[i] <- (bc.low <= pi & pi <= bc.high)
    bc.cover0[i] <- (bc.low <= 0)
    bc.cover.direct[i] <- (bc.low <= direct & direct <= bc.high)
    REbc[i] <- (bc.high - bc.low) / (naive.high - naive.low)
  }
  close(pb)

  #############################################################
  # COMPUTE SUMMARY STATISTICS
  #############################################################

  # Average bias
  bias.naive.cm <- mean(naive.pred - pi)
  bias.bias.correct <- mean(bc.pred - pi)

  # Root Mean Squared Error
  rmse.naive.cm <- mean((naive.pred - pi)^2)
  rmse.bias.correct <- mean((bc.pred - pi)^2)

  # Coverage rates
  coverage.naive.cm <- mean(naive.cover, na.rm = TRUE) * 100
  coverage.bias.correct <- mean(bc.cover, na.rm = TRUE) * 100
  coverage.zero <- mean(bc.cover0, na.rm = TRUE) * 100
  coverage.direct <- mean(bc.cover.direct, na.rm = TRUE) * 100

  # Sort estimates for visualization purposes
  sorted_idx <- order(bc.pred)
  bc.pred.sort <- bc.pred[sorted_idx]
  bc.low.save.sort <- bc.low.save[sorted_idx]
  bc.high.save.sort <- bc.high.save[sorted_idx]

  #########
  # RETURN OUTPUT
  #########
  outvec <- c(bias.naive.cm, bias.bias.correct,
              rmse.naive.cm, rmse.bias.correct,
              coverage.naive.cm, coverage.bias.correct,
              coverage.zero, coverage.direct)
  outvec <- sapply(outvec, roundy)
  names(outvec) <- c(
    "Average Bias (naive)", "Average Bias (bias-corrected)",
    "RMSE (naive)", "RMSE (bias-corrected)",
    "% Coverage (naive)", "% Coverage (bias-corrected)",
    "% Coverage (bias-corrected) of 0",
    "% Coverage (bias-corrected) of Direct Estimate"
  )

  out <- list(
    Results = outvec,
    BiasCorrectEst = bc.pred.sort,
    BiasCorrectLow = bc.low.save.sort,
    BiasCorrectHigh = bc.high.save.sort,
    EstimatedBias = est.bias,
    RelativeLengthCI = REbc
  )

  return(out)
}
