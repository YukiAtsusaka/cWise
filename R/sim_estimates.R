#' Simulate and Plot Panel C of Figure C7
#'
#' Runs Monte Carlo simulations using the bias-corrected crosswise model and
#' creates a caterpillar plot showing sorted point estimates with bootstrap
#' confidence intervals, replicating Panel C of Figure C7 in Atsusaka and
#' Stevenson (2021).
#'
#' @param N.sim Integer. Number of Monte Carlo simulations. Default is 100.
#' @param sample Integer. Sample size per simulation.
#' @param pi Numeric. True prevalence rate of the sensitive attribute.
#' @param p Numeric. Randomization probability for the sensitive question.
#' @param p.prime Numeric. Randomization probability for the anchor question.
#' @param gamma Numeric. Proportion of attentive respondents (between 0 and 1).
#' @param direct Numeric. Direct questioning estimate for comparison.
#' @param txcol Character. Color for annotation text. Default: \code{"dimgray"}.
#' @param sim.results Optional list. Pre-computed output from \code{\link{sim.cwdata}}.
#'   If \code{NULL} (default), the simulation is run internally.
#'
#' @return Invisibly returns the simulation results list from \code{\link{sim.cwdata}},
#'   containing \code{BiasCorrectEst}, \code{BiasCorrectLow},
#'   \code{BiasCorrectHigh}, and summary \code{Results}.
#'
#' @details
#' The plot displays:
#' \itemize{
#'   \item Sorted bias-corrected point estimates as filled circles
#'   \item Bootstrap 95\% confidence intervals as vertical line segments
#'   \item A horizontal reference line at 0
#'   \item A horizontal reference line at the true prevalence \code{pi} (red)
#' }
#' Estimates are sorted in ascending order, creating a characteristic
#' "fan" shape that reveals the distribution of estimates across simulations.
#'
#' Text annotation positions are calibrated to \code{N.sim = 100} and scale
#' proportionally for other values.
#'
#' If \code{sim.results} is supplied, all simulation parameters (\code{N.sim},
#' \code{sample}, \code{p}, \code{p.prime}, \code{gamma}, \code{direct}) are
#' still used for the annotations, but no new simulation is run.
#'
#' @examples
#' \dontrun{
#' # Replicate Panel C of Figure C7
#' sim.estimates(
#'   N.sim   = 100,
#'   sample  = 1000,
#'   pi      = 0.1,
#'   p       = 0.1,
#'   p.prime = 0.1,
#'   gamma   = 0.8,
#'   direct  = 0.1
#' )
#'
#' # Re-use pre-computed simulation results
#' res <- sim.cwdata(N.sim = 100, sample = 1000, pi = 0.1,
#'                   p = 0.1, p.prime = 0.1, gamma = 0.8, direct = 0.1)
#' sim.estimates(sample = 1000, pi = 0.1, p = 0.1, p.prime = 0.1,
#'               gamma = 0.8, direct = 0.1, sim.results = res)
#' }
#'
#' @seealso \code{\link{sim.cwdata}} for the underlying simulation function
#'
#' @references
#' Atsusaka and Stevenson (2021). Figure C7, Panel C.
#'
#' @export
sim.estimates <- function(N.sim = 100, sample, pi, p, p.prime, gamma, direct,
                          txcol = "dimgray", sim.results = NULL) {

  # Run simulation if pre-computed results are not supplied
  if (is.null(sim.results)) {
    sim.results <- sim.cwdata(
      N.sim   = N.sim,
      sample  = sample,
      pi      = pi,
      p       = p,
      p.prime = p.prime,
      gamma   = gamma,
      direct  = direct
    )
  }

  est    <- sim.results$BiasCorrectEst
  low    <- sim.results$BiasCorrectLow
  high   <- sim.results$BiasCorrectHigh
  n_sims <- length(est)

  # Scale factor so text annotations work for any N.sim (calibrated for N.sim = 100)
  sf <- n_sims / 100

  # Caterpillar plot: sorted estimates with bootstrap CIs
  plot(est, type = "n", ylim = c(0, max(high)), xlab = "", las = 1, ylab = "")
  points(est, pch = 16, col = scales::alpha("dimgray", 1))
  arrows(
    x0 = seq_along(est), x1 = seq_along(est),
    y0 = low,            y1 = high,
    length = 0, col = scales::alpha("gray60", 0.9)
  )
  abline(h = 0,  col = "dimgray",    lwd = 1.5)
  abline(h = pi, col = "firebrick4", lty = 1, lwd = 1.5)

  # legend(
  #   x = 0, y = 0.42,
  #   legend = "Quantity of interest",
  #   lty = 1, lwd = 1.5, col = "firebrick4",
  #   box.lty = 0, cex = 1.5
  # )

  # # Parameter annotations — x positions scale with n_sims
  # text(x = 18   * sf, y = 0.33, labels = bquote("n"          == .(sample)),  cex = 1.5, col = txcol)
  # text(x = 15   * sf, y = 0.30, labels = bquote(pi           == .(pi)),      cex = 1.5, col = txcol)
  # text(x = 12.5 * sf, y = 0.27, labels = bquote(paste(pi, "'") == 0),        cex = 1.5, col = txcol)
  # text(x = 15   * sf, y = 0.24, labels = bquote("p"          == .(p)),       cex = 1.5, col = txcol)
  # text(x = 16   * sf, y = 0.21, labels = bquote("p'"         == .(p.prime)), cex = 1.5, col = txcol)
  # text(x = 14.5 * sf, y = 0.18, labels = bquote(gamma        == .(gamma)),   cex = 1.5, col = txcol)

  mtext(side = 1, "Simulation Index",                  line = 2)
  # mtext(side = 2, "Prevalence of Sensitive Attribute", line = 2)
  title("Simulated Bias-Corrected Estimates",
        adj = 0, cex.main = 1, font = 2, line = 0.5)
  box()

  invisible(sim.results)
}
