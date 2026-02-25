#' Determine Sample Size for Desired Coverage Properties
#'
#' Tests multiple sample sizes to determine which achieves desired confidence
#' interval coverage properties. Specifically, it calculates what percentage of
#' 95\% confidence intervals (1) exclude zero and (2) include the direct estimate.
#' This helps researchers choose a sample size that provides sufficient precision.
#'
#' @param N.sim Integer. Number of Monte Carlo simulations per sample size.
#'   Default is 50. Larger values provide more stable estimates but increase
#'   computation time.
#' @param pi Numeric. True prevalence rate of the sensitive attribute (between 0 and 1).
#' @param p Numeric. Probability for the randomization item in sensitive question.
#' @param p.prime Numeric. Probability for the anchor question (non-sensitive).
#' @param gamma Numeric. Proportion of attentive respondents (between 0 and 1).
#' @param direct Numeric. Direct questioning estimate for comparison purposes.
#'
#' @return A data frame with three columns:
#' \describe{
#'   \item{SampleSize}{Sample sizes tested: 100, 500, 1000, 1500, 2000, 2500, 3000}
#'   \item{CoverageZero}{Percentage of 95\% CIs that include zero.
#'                       Lower values indicate better precision (CIs exclude zero).}
#'   \item{CoverageDirect}{Percentage of 95\% CIs that include the direct estimate.
#'                         Values near 95\% suggest good agreement with direct questioning.}
#' }
#'
#' @details
#' This function is useful for planning studies where researchers want to:
#' \enumerate{
#'   \item Distinguish the estimated prevalence from zero with high confidence
#'   \item Obtain narrow confidence intervals for precise estimation
#'   \item Compare crosswise estimates with direct questioning estimates
#' }
#'
#' For each sample size, the function:
#' \itemize{
#'   \item Simulates N.sim datasets using the crosswise model
#'   \item Computes bias-corrected estimates with bootstrap 95\% CIs
#'   \item Calculates what percentage of CIs contain zero
#'   \item Calculates what percentage of CIs contain the direct estimate
#' }
#'
#' A progress bar displays simulation progress.
#'
#' @note
#' \itemize{
#'   \item Low CoverageZero values indicate CIs that reliably exclude zero
#'         (good for establishing that prevalence is non-zero)
#'   \item CoverageDirect near 95\% suggests consistency between crosswise and
#'         direct questioning approaches
#'   \item This function uses 500 bootstrap iterations per simulation
#' }
#'
#' @examples
#' # Find sample size needed to reliably exclude zero
#' \dontrun{
#' result <- sim.power.N(
#'   N.sim = 50,
#'   pi = 0.1,
#'   p = 0.1,
#'   p.prime = 0.1,
#'   gamma = 0.8,
#'   direct = 0.05
#' )
#'
#' print(result)
#'
#' # Visualize results
#' plot(result$SampleSize, result$CoverageZero,
#'      type = "b", xlab = "Sample Size",
#'      ylab = "% of CIs Including Zero",
#'      main = "Precision vs Sample Size")
#' }
#'
#' @seealso \code{\link{sim.power}} for the underlying simulation function
#'
#' @references
#' Atsusaka and Stevenson (2021). Appendix C5: Sample Size Determination
#' and Parameter Selection.
#'
#' @export
sim.power.N <- function(N.sim = 50, pi, p, p.prime, gamma, direct) {

  # Helper function for rounding
  roundy <- function(x) { round(x, digits = 3) }

  p2 <- p.prime
  Nvec <- c(100, 500, 1000, 1500, 2000, 2500, 3000)
  pb <- txtProgressBar(min = 1, max = length(Nvec), initial = 0, style = 3)

  coverage.zero.save <- numeric(length(Nvec))
  coverage.direct.save <- numeric(length(Nvec))

  # Loop through each sample size
  for (n in 1:length(Nvec)) {

    setTxtProgressBar(pb, n)
    N <- Nvec[n]

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
    # Repeat simulation N.sim times with current sample size
    #############################################################
    for (i in 1:N.sim) {

      # Generate attentive/inattentive status
      Attentive <- rbinom(n = N, size = 1, prob = gamma)

      # SENSITIVE QUESTION OF INTEREST (Crosswise format)
      StatementA <- rbinom(n = N, size = 1, prob = pi)           # Sensitive item
      StatementB <- rbinom(n = N, size = 1, prob = p)            # Randomization item
      True.res <- ifelse(StatementA != StatementB, 0, 1)         # True answer

      Obs.res <- rep(NA, N)
      Obs.res[Attentive == 1] <- True.res[Attentive == 1]        # Attentive Rs
      Obs.res[Attentive == 0] <- rbinom(n = sum(Attentive == 0),
                                        size = 1, prob = 0.5)     # Inattentive Rs

      # NON-SENSITIVE AUXILIARY QUESTION (Anchor question)
      StatementC <- 0                                             # Known prevalence
      StatementD <- rbinom(n = N, size = 1, prob = p2)            # Randomization
      True.res.prime <- ifelse(StatementC != StatementD, 0, 1)

      Obs.res.prime <- rep(NA, N)
      Obs.res.prime[Attentive == 1] <- True.res.prime[Attentive == 1]
      Obs.res.prime[Attentive == 0] <- rbinom(n = sum(Attentive == 0),
                                              size = 1, prob = 0.5)

      dat <- data.frame(Obs.res, True.res, Attentive, Obs.res.prime, True.res.prime)


      # (2) NAIVE CROSSWISE MODEL ESTIMATOR
      lambda.hat <- mean(Obs.res)
      pi.hat.naive <- (lambda.hat + p - 1) / (2*p - 1)

      pi.hat.naive.var <- (pi.hat.naive * (1 - pi.hat.naive)) / (N - 1) +
                          (p * (1 - p)) / ((N - 1) * ((2*p - 1)^2))
      pi.hat.naive.sd <- sqrt(pi.hat.naive.var)

      naive.low  <- max(0, pi.hat.naive - 1.96 * pi.hat.naive.sd)
      naive.high <- min(1, pi.hat.naive + 1.96 * pi.hat.naive.sd)


      # (3) BIAS CORRECTED CROSSWISE MODEL ESTIMATOR
      gamma.hat <- (mean(Obs.res.prime) - 0.5) / (0.5 - p2)
      Bias.hat <- (1/2) * ((lambda.hat - 0.5) / (p - 0.5)) -
                  (1 / (2 * gamma.hat)) * ((lambda.hat - 0.5) / (p - 0.5))

      pi.hat.bc <- pi.hat.naive - Bias.hat
      pi.hat.bc <- max(0, min(1, pi.hat.bc))

      # BOOTSTRAPPING (500 iterations)
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
        bs[b] <- max(0, min(1, bs[b]))
      }

      bc.low  <- max(0, quantile(bs, prob = 0.025, na.rm = TRUE))
      bc.high <- min(1, quantile(bs, prob = 0.975, na.rm = TRUE))


      # STORE RESULTS
      est.bias[i] <- Bias.hat
      naive.pred[i] <- pi.hat.naive
      bc.pred[i] <- pi.hat.bc
      bc.high.save[i] <- bc.high
      bc.low.save[i] <- bc.low
      naive.cover[i] <- (naive.low <= pi & pi <= naive.high)
      bc.cover[i] <- (bc.low <= pi & pi <= bc.high)
      bc.cover0[i] <- (bc.low <= 0)                              # Does CI include 0?
      bc.cover.direct[i] <- (bc.low <= direct & direct <= bc.high)  # Does CI include direct estimate?
      REbc[i] <- (bc.high - bc.low) / (naive.high - naive.low)
    }

    #############################################################
    # COMPUTE COVERAGE STATISTICS FOR THIS SAMPLE SIZE
    #############################################################
    coverage.zero <- mean(bc.cover0, na.rm = TRUE) * 100
    coverage.direct <- mean(bc.cover.direct, na.rm = TRUE) * 100

    coverage.zero.save[n] <- coverage.zero
    coverage.direct.save[n] <- coverage.direct
  }
  close(pb)

  #########
  # RETURN OUTPUT
  #########
  out <- data.frame(
    SampleSize = Nvec,
    CoverageZero = coverage.zero.save,
    CoverageDirect = coverage.direct.save
  )

  return(out)
}
