#' @title cmBound
#'
#' @description \code{cmBound} is used to apply a sensitivity analysis and
#' visualize the sensitivity bounds for naive crosswise estimates. The
#' sensitivity analysis assumes that inattentive respondents choose either
#' response option with equal probability.
#'
#' @param lambda.hat a value for the observed crosswise proportion: Prop(TRUE-TRUE or FALSE-FALSE).
#' @param p an auxiliary probability for the crosswise question. It must be a
#' single finite value strictly between 0 and 1, excluding 0.5.
#' @param N a value for the number of survey respondents in the crosswise (and anchor) question.
#' @param dq a value for a point estimate from direct questioning (if available).
#' @param N.dq a value for the number of survey respondents in direct questioning (if available).
#'
#' @return A ggplot object visualizing the result of sensitivity analysis.
#' @seealso [bc_est()] for bias-corrected prevalence estimation.
#' @references Atsusaka, Y. and Stevenson, R. T. (2023). The crosswise model
#' for sensitive survey questions. \doi{10.1017/pan.2021.43}.
#' @examples
#' p <- cmBound(lambda.hat=0.6385, p=0.25, N=310, dq=0.073, N.dq=310)
#' p
#'
#' p <- p + ggplot2::ggtitle("Sensitivity Analysis") +
#'          ggplot2::theme(plot.title = ggplot2::element_text(size=20, face="bold"))
#' @export
#' @import ggplot2


cmBound = function(lambda.hat, p, N, dq=NULL, N.dq=NULL){

  if (length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1 || p == 0.5) {
    stop(
      "`p` must be a single finite probability strictly between 0 and 1, excluding 0.5.",
      call. = FALSE
    )
  }

  pi.hat.naive = (lambda.hat+p-1)/(2*p-1)
  pi.hat.naive.var = (pi.hat.naive*(1-pi.hat.naive))/(N-1) + (p*(1-p))/((N-1)*((2*p-1)^2))
  pi.hat.naive.sd = sqrt(pi.hat.naive.var)

  gamma.hat = rev(seq(from=0.5, to=1, by=0.01)) # Vector of attentitive rate
  Bias.hat = (1/2)*((lambda.hat-0.5)/(p-0.5)) - (1/(2*gamma.hat))*((lambda.hat-0.5)/(p-0.5))

  pi.hat.bc = rep(pi.hat.naive, length(gamma.hat)) - Bias.hat # Bias Correction
  pi.hat.bc = ifelse(pi.hat.bc < 1, pi.hat.bc, 1) # Logical bound
  pi.hat.bc = ifelse(pi.hat.bc > 0, pi.hat.bc, 0) # Logical bound
  pi.hat.bc.var = (1/(gamma.hat^2))*
    ((pi.hat.bc*(1-pi.hat.bc))/(N-1) +
       (p*(1-p))/((N-1)*((2*p-1)^2)))
  pi.hat.bc.sd = sqrt(pi.hat.bc.var)

  mean = pi.hat.bc
  low  = pi.hat.bc-1.96*pi.hat.bc.sd; low  = ifelse(low > 0, low, 0)
  high = pi.hat.bc+1.96*pi.hat.bc.sd; high = ifelse(high < 1, high, 1)

  perinat = (1-gamma.hat)*100
  ggdata <- data.frame(mean, low, high, perinat)

  plot_object <- ggplot(ggdata, aes(x=perinat, y=mean)) +
       geom_line(col="firebrick", linewidth=1.2) +
       geom_ribbon(aes(ymin=low, ymax=high),linetype=2, alpha=0.5)+
       geom_point(x=0, y=mean[1], shape=4, size=6, color="navy")+
       xlab("Inattentive Respondents (%)")+
       ylab("Estimated Proportion")+ ylim(0,0.5)+
       theme_bw()+
       theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             plot.title = element_text(size = 10, face = "bold"),
             axis.text.x = element_text(size = 18),
             axis.title.x = element_text(size = 18),
             axis.text.y = element_text(size = 15),
             axis.title.y = element_text(size = 15))

  if (!is.null(dq)) {
    if (length(dq) != 1L || !is.finite(dq) || dq < 0 || dq > 1 ||
        length(N.dq) != 1L || !is.finite(N.dq) || N.dq <= 1) {
      stop("`dq` must be in [0, 1] and `N.dq` must be greater than 1.",
           call. = FALSE)
    }
    dq.var <- dq * (1 - dq) / (N.dq - 1)
    dq.sd <- sqrt(dq.var)
    dq.upper <- min(dq + 1.96 * dq.sd, 1)
    dq.lower <- max(dq - 1.96 * dq.sd, 0)
    plot_object <- plot_object +
      geom_hline(yintercept = dq, col = "dimgray", linewidth = 0.9) +
      geom_hline(yintercept = dq.upper, linetype = "dashed", col = "dimgray", linewidth = 0.9) +
      geom_hline(yintercept = dq.lower, linetype = "dashed", col = "dimgray", linewidth = 0.9)
  }

  return(plot_object)
}
