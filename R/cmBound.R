#' @title cmBound
#'
#' @description Visualize the sensitivity bounds for naive crosswise estimates
#'
#' @param lambda.hat observed crosswise proportion: Prop(TRUE-TRUE or FALSE-FALSE)
#' @param p Known proportion for the non-sensitive statement
#' @param N Number of survey respondens in direct questioning (if available)
#' @param dq Point estimate from direct questioning (if available)
#'
#' @return A ggplot object.
#' @examples
#' p <- cmBound(lambda.hat=0.6385, p=0.25, N=310, dq=0.073)
#' p
#'
#' p <- p + ggtitle("Sensitivity Analysis") +
#'          theme(plot.title = element_text(size=20, face="bold"))
#' @export
#' @importFrom ggplot2


cmBound = function(lambda.hat, p, N, dq=NULL, N.dq=NULL){

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

  dq.var = dq*(1-dq)/(N.dq-1)
  dq.sd = sqrt(dq.var)
  dq.upper = min(dq+1.96*dq.sd,1)
  dq.lower = max(dq-1.96*dq.sd,0)

  perinat = (1-gamma.hat)*100
  ggdata <- as.data.frame(cbind(mean, low, high, perinat, dq))

  p <- ggplot(ggdata, aes(x=perinat, y=mean)) +
       geom_hline(yintercept=dq, col="dimgray", size=1.2) +
       geom_hline(yintercept=dq.upper, linetype="dashed", col="dimgray", size=1.2) +
       geom_hline(yintercept=dq.lower, linetype="dashed", col="dimgray", size=1.2) +
       geom_line(col="firebrick", size=2) +
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
  return(p)
}
