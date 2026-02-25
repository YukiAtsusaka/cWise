#' Generate Power Curves for Sample Size Determination
#'
#' Computes statistical power across a range of sample sizes for a one-sided
#' hypothesis test of H0: π ≤ π0 vs H1: π > π0. This function is used to create
#' power curves (Panel A in Figure C7) to help researchers determine the sample
#' size needed to achieve desired statistical power.
#'
#' @param N.sim Integer. Number of Monte Carlo simulations per sample size.
#'   Larger values provide more stable power estimates but increase computation time.
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
#' @return A list containing:
#' \describe{
#'   \item{simulated_power}{Numeric vector of statistical power values corresponding
#'                          to each sample size tested}
#'   \item{sample_size}{Numeric vector of sample sizes tested: c(1, 500, 1000, 1500, 2000, 2500)}
#' }
#'
#' @details
#' The function implements the power calculation based on the Wald test:
#' \deqn{Power = \beta = 1 - \Phi\left(\frac{\pi_0 - \pi_1 - c\tilde{\sigma}_0}{\tilde{\sigma}_1}\right)}
#' where:
#' \itemize{
#'   \item \eqn{\Phi} is the cumulative distribution function of the standard normal
#'   \item \eqn{c = \Phi^{-1}(1-\alpha)} is the critical value for significance level \eqn{\alpha}
#'   \item \eqn{\tilde{\sigma}_0} and \eqn{\tilde{\sigma}_1} are simulated standard errors
#'         under H0 and H1 respectively
#' }
#'
#' The function tests sample sizes of 1, 500, 1000, 1500, 2000, and 2500. For each
#' sample size, it runs sim.power() twice (once under H0 and once under H1) to
#' estimate the sampling distributions.
#'
#' @note This function can take considerable time to run (potentially hours) depending
#' on N.sim. For quick exploration, use smaller N.sim values (e.g., 100-500). For
#' publication-quality results, use N.sim >= 2000.
#'
#' @examples
#' # Generate power curve for p = 0.1
#' \dontrun{
#' curve_result <- sim.curve(
#'   N.sim = 500,
#'   pi.null = 0,
#'   pi.alt = 0.1,
#'   p = 0.1,
#'   p.prime = 0.1,
#'   gamma = 0.8,
#'   direct = 0.02
#' )
#'
#' # Plot the power curve
#' plot(curve_result$sample_size, curve_result$simulated_power,
#'      type = "l", xlab = "Sample Size", ylab = "Power",
#'      main = "Statistical Power vs Sample Size")
#' abline(h = 0.8, lty = 2, col = "red")  # 80% power line
#' }
#'
#' # To find sample size for 80% power, look where curve crosses 0.8
#'
#' @seealso \code{\link{sim.power}} for the underlying simulation function
#'
#' @references
#' Ulrich et al. (2012). Using the weighted likelihood ratio test for composite
#' hypotheses in diagnostic studies. Statistics in Medicine.
#'
#' Atsusaka and Stevenson (2021). Appendix C5: Sample Size Determination
#' and Parameter Selection.
#'
#' @export
sim.curve <- function(N.sim, pi.null, pi.alt, p, p.prime, gamma, direct) {

  # Sample sizes to test
  n_vec <- c(1, seq(from = 500, to = 2500, by = 500))
  out_n <- numeric(length(n_vec))
  out_power <- numeric(length(n_vec))

  # Significance level
  alpha <- 0.05
  mu0 <- pi.null
  mu1 <- pi.alt

  # Loop through each sample size
  for (i in 1:length(n_vec)) {

    n <- n_vec[i]
    cat(sprintf("Now, working on the sample size of %d\n", n))

    # Simulate under H0 (null hypothesis)
    H0 <- sim.power(
      N.sim = N.sim,
      sample = n,
      pi = pi.null,
      p = p,
      p.prime = p.prime,
      gamma = gamma,
      direct = direct
    )

    # Simulate under H1 (alternative hypothesis)
    H1 <- sim.power(
      N.sim = N.sim,
      sample = n,
      pi = pi.alt,
      p = p,
      p.prime = p.prime,
      gamma = gamma,
      direct = direct
    )

    # CALCULATE STATISTICAL POWER
    # Based on formula: Power = 1 - Φ((π0 - π1 - c*σ0) / σ1)
    # where Φ is the standard normal CDF

    # Simulated standard errors
    sim.sigma0 <- sd(H0$BiasCorrectEst)  # SD under H0
    sim.sigma1 <- sd(H1$BiasCorrectEst)  # SD under H1

    # Critical value calculation
    c_one_minus_alpha <- quantile(H0$BiasCorrectEst, 1 - alpha)
    z_one_minus_alpha <- (c_one_minus_alpha - mu0) / sim.sigma0
    z_alpha <- -1 * z_one_minus_alpha

    # Power calculation
    A <- (mu1 - mu0 + z_alpha * sim.sigma0) / sim.sigma1
    power <- pnorm(A, mean = 0, sd = 1)

    # Store results
    out_power[i] <- power
    out_n[i] <- n
  }

  # Return results
  out <- list(
    simulated_power = out_power,
    sample_size = out_n
  )

  return(out)
}
