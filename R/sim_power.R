#' Compute Statistical Power for a Fixed Sample Size
#'
#' Computes statistical power for a one-sided hypothesis test of
#' H0: π ≤ π0 vs H1: π > π0 at a given fixed sample size. This is useful
#' when researchers already know their sample size and want to assess the
#' power they can expect from their study.
#'
#' @param N.sim Integer. Number of Monte Carlo simulations.
#'   Larger values provide more stable power estimates but increase computation time.
#' @param sample Integer. The fixed sample size (number of respondents) to evaluate.
#' @param pi.null Numeric. Prevalence rate under the null hypothesis (π0).
#'   Can be 0 or a value from direct questioning.
#' @param pi.alt Numeric. True prevalence rate under the alternative hypothesis (π1).
#'   Must be greater than pi.null.
#' @param p Numeric. Probability for the randomization item in the sensitive question.
#'   Values between 0.1 and 0.3 are typical.
#' @param p.prime Numeric. Probability for the anchor question (non-sensitive).
#' @param gamma Numeric. Proportion of attentive respondents (between 0 and 1).
#'   For example, 0.8 means 80\% of respondents are attentive.
#' @param direct Numeric. Direct questioning estimate for comparison purposes.
#'
#' @return A numeric scalar representing the estimated statistical power
#'   (probability of correctly rejecting H0 when H1 is true) at the given sample size.
#'
#' @details
#' The function implements the power calculation based on the Wald test:
#' \deqn{Power = \Phi\left(\frac{\pi_1 - \pi_0 + z_\alpha \tilde{\sigma}_0}{\tilde{\sigma}_1}\right)}
#' where:
#' \itemize{
#'   \item \eqn{\Phi} is the cumulative distribution function of the standard normal
#'   \item \eqn{z_\alpha} is derived from the simulated null distribution
#'   \item \eqn{\tilde{\sigma}_0} and \eqn{\tilde{\sigma}_1} are simulated standard errors
#'         under H0 and H1 respectively
#' }
#'
#' The function runs \code{sim.cwdata()} twice — once under H0 using \code{pi.null}
#' and once under H1 using \code{pi.alt} — to estimate the sampling distributions
#' at the specified sample size.
#'
#' @note This function can take considerable time to run depending on N.sim and
#' sample size. For quick exploration, use smaller N.sim values (e.g., 100–500).
#' For publication-quality results, use N.sim >= 2000.
#'
#' @examples
#' # Compute power at a fixed sample size of 1000
#' \dontrun{
#' pwr <- sim.power(
#'   N.sim  = 500,
#'   sample = 1000,
#'   pi.null = 0,
#'   pi.alt  = 0.1,
#'   p       = 0.1,
#'   p.prime = 0.1,
#'   gamma   = 0.8,
#'   direct  = 0.02
#' )
#' cat(sprintf("Estimated power: %.3f\n", pwr))
#' }
#'
#' @seealso \code{\link{sim.cwdata}} for the underlying simulation function
#'
#' @references
#' Atsusaka and Stevenson (2021). Appendix C5: Sample Size Determination
#' and Parameter Selection.
#'
#' @export
sim.power <- function(N.sim, sample, pi.null, pi.alt, p, p.prime, gamma, direct) {

  alpha <- 0.05
  mu0   <- pi.null
  mu1   <- pi.alt

  cat(sprintf("Simulating under H0 (pi = %.3f, n = %d)...\n", mu0, sample))

  # Simulate under H0 (null hypothesis)
  H0 <- sim.cwdata(
    N.sim   = N.sim,
    sample  = sample,
    pi      = pi.null,
    p       = p,
    p.prime = p.prime,
    gamma   = gamma,
    direct  = direct
  )

  cat(sprintf("Simulating under H1 (pi = %.3f, n = %d)...\n", mu1, sample))

  # Simulate under H1 (alternative hypothesis)
  H1 <- sim.cwdata(
    N.sim   = N.sim,
    sample  = sample,
    pi      = pi.alt,
    p       = p,
    p.prime = p.prime,
    gamma   = gamma,
    direct  = direct
  )

  # Simulated standard errors
  sim.sigma0 <- sd(H0$BiasCorrectEst)
  sim.sigma1 <- sd(H1$BiasCorrectEst)

  # Critical value derived from the null distribution
  c_one_minus_alpha  <- quantile(H0$BiasCorrectEst, 1 - alpha)
  z_one_minus_alpha  <- (c_one_minus_alpha - mu0) / sim.sigma0
  z_alpha            <- -1 * z_one_minus_alpha

  # Power: P(reject H0 | H1 true)
  A     <- (mu1 - mu0 + z_alpha * sim.sigma0) / sim.sigma1
  power <- pnorm(A, mean = 0, sd = 1)

  return(power)
}
