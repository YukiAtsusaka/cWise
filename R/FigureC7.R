#############################################################################
# FigureC7.R
# Aim: to replicate Figure C7 in Atsusaka and Stevenson (2021)
# Run time: about 2.5 HOURS
#############################################################################

rm(list=ls())
start <- Sys.time()
library(tidyverse)
source("R/FigureC7_function.R") # Read Three Functions

############################################################
# FOR PANEL A
############################################################

curve_list1 <- list()
#curve_list2 <- list()
p_vec <- c(0.1, 0.2, 0.3) # Non-sensitive Probability

# Setting Gamma=0.8
for(k in 1:3){
curve_list1[[k]] <- sim.curve(N.sim=2000, pi.null=0, pi.alt=0.1,
                             p=p_vec[k], # VARYING THIS PARAMETER
                             p.prime=0.1, gamma=0.8, direct=0.02)
}


############################################################
# FOR PANEL B
############################################################

N.sim <- 2000
set.seed(1235)
# SIMULATING POINT ESTIMATES UNDER H0 (NULL) AND H1 (ALTERNATIVE HYPOTHESIS)
H0 <- sim.power(N.sim=2000, sample=400, pi=0, p=0.1, p.prime=0.1, gamma=0.8, direct=0.1)
H1 <- sim.power(N.sim=2000, sample=400, pi=0.1, p=0.1, p.prime=0.1, gamma=0.8, direct=0.1)

c_one_minus_alpha <- quantile(H0$BiasCorrectEst, 0.95) # Critical Value under the Normal Approximation
indx <- ifelse(H1$BiasCorrectEst>c_one_minus_alpha,T,F)
digit <- 1:N.sim

y_H1 <- ecdf(H1$BiasCorrectEst)
y_H1 <- y_H1(H1$BiasCorrectEst)
y_axisH1 <- y_H1[min(digit[indx==T])]

y_H0 <- ecdf(H0$BiasCorrectEst)
y_H0 <- y_H0(H0$BiasCorrectEst)

###############
# TYPE I ERROR
###############
indx0 <- ifelse(H0$BiasCorrectEst>c_one_minus_alpha,T,F)
est_max <-  max(H0$BiasCorrectEst)
y_H0 <-  y_H0[y_H0 > 0.95]
TypeI_x <- c(c_one_minus_alpha,
             H0$BiasCorrectEst[min(digit[indx0==T]):N.sim],
             est_max)
TypeI_y <- c(0.95, y_H0, 0.95)


####################
# STATISTICAL POWER
####################

est_max <-  max(H1$BiasCorrectEst)
y_H1_power <-  y_H1[y_H1 >= y_axisH1]
Power_x <- c(c_one_minus_alpha,
             H1$BiasCorrectEst[min(digit[indx==T]):N.sim],
             est_max)
Power_y <- c(y_axisH1, y_H1_power, y_axisH1)


#################
# TYPE II ERROR
#################
y_H1_typeII <- y_H1[y_H1 < y_axisH1]
TypeII_a_x <- c(0,
                H1$BiasCorrectEst[1:max(digit[indx==F])],
                c_one_minus_alpha)
TypeII_a_y <- c(0, y_H1_typeII, 0)

# RECTANGULAR PART
TypeII_b_x <- c(c_one_minus_alpha, c_one_minus_alpha, est_max, est_max)
TypeII_b_y <- c(0, y_axisH1, y_axisH1, 0)


############################################################
# FOR PANEL C
############################################################

test <- sim.power(N.sim=100,
                  sample=1000,
                  pi=0.1,
                  p=0.1,
                  p.prime=0.1,
                  gamma=0.8,
                  direct=0.1)

test$Results

###################################################
# VISUALIZE TYPE I ERROR, POWER, AND TYPE II ERROR
###################################################

library(rlist)
list.save(curve_list1, 'list1.rdata')

load('list1.rdata')
curve_list1 <- x

#############################################################################
# CREATE A SINGLE GRAPH
#############################################################################

pdf("FigureC7.pdf", width=9, height=3.2)
par(mfrow=c(1,3), mar=c(4,4,2,1), oma=c(0,0,1,0)+0.01, mgp=c(3,0.6,0))

############################################################
# PANEL A
############################################################
txcol = "dimgray"
# VISUALIZE THE POWER CURVE FOR VARIOUS CONDITIONS
        plot(0, 0, xlim=c(0,2500), ylim=c(0,1), las=1, type="n", yaxt="n", xaxt="n", xlab="", ylab="")
        axis(side=2, at=seq(from=0,to=1,by=0.1),las=1)
        axis(side=1, at=seq(from=0,to=2500,by=500))
        lty_vec <- c(1,2,3)
        # DRAWING LINES
        for(k in 1:3){
        power <- curve_list1[[k]]$simulated_power
        size <- curve_list1[[k]]$sample_size
        lines(size, power, lwd=2, lty=lty_vec[k], col="firebrick4")
        }

legend(x=250,y=0.35, legend=c("p = 0.1", "p = 0.2", "p = 0.3"), cex=1.6,
       lty=c(1,2,3), lwd=2, col="firebrick4", box.lty=0)

text(x=2000, y=0.8, labels= bquote(pi[0] == 0), cex=1.5, col=txcol)
text(x=2000+90, y=0.7, labels= bquote(pi[1] == 0.1), cex=1.5, col=txcol)
text(x=2000, y=0.6, labels= bquote(paste(pi,"'") == 0), cex=1.5, col=txcol)
text(x=2000+80, y=0.5, labels= bquote("p'" == 0.1), cex=1.5, col=txcol)
text(x=2000+50, y=0.4, labels= bquote(gamma == 0.8), cex=1.5, col=txcol)
text(x=2000+100, y=0.3, labels= bquote(alpha == 0.05), cex=1.5, col=txcol)

# AXIS TITLES
mtext(side=1, "Sample Size", line=2)
mtext(side=2, "Power", line=2)
title("A", adj=0, cex.main=2, font=2, line=0.5)
box()

############################################################
# PANEL B
############################################################

plot(ecdf(H1$BiasCorrectEst), col="white", cex=0, main="",
     ylab="", xlab="", las=1, xlim=c(0,0.22))
abline(v=c_one_minus_alpha, col="gray50", lwd=1, lty=2) # Vertical Line
polygon(Power_x, Power_y,col=gray(0.9), border=NA)                     # POWER
polygon(TypeI_x,TypeI_y,col=alpha("navy",0.9), border=NA)              # TYPE I ERROR
polygon(TypeII_a_x,TypeII_a_y,col=alpha("firebrick4",0.9), border=NA)  # TYPE II ERROR (a)
polygon(TypeII_b_x,TypeII_b_y,col=alpha("firebrick4",0.9), border=NA)  # TYPE II ERROR (b)
lines(ecdf(H1$BiasCorrectEst), col="firebrick4", cex=0, lwd=1.5)
lines(ecdf(H0$BiasCorrectEst), col="dimgray", cex=0, lwd=1.5)

text(x=0.0625, y=0.85, labels="Type I Error", font=2, cex=1.8, col="navy")
text(x=0.15, y=0.3, labels="Power", font=2, cex=2, col="dimgray")
text(x=0.125, y=0.08, labels="Type II Error", font=2, cex=1.8, col="white")
text(x=0.015,y=0.45, labels=expression(c[1-alpha]), font=1.8, col="dimgray", cex=1.5)
arrows(x0=0.035,x1=c_one_minus_alpha-0.002,y0=0.45,y1=0.45,length=0.08, col="dimgray")
text(x=0.03, y=0.37, labels= bquote(alpha == 0.05), cex=1.3, col=txcol, font=2)

text(x=0.17+0.006, y=0.84, labels= bquote("n" == 400), cex=1.5, col=txcol)
text(x=0.17+0.005, y=0.77, labels= bquote("p" == 0.1), cex=1.5, col=txcol)
text(x=0.17, y=0.70, labels= bquote(pi[0] == 0), cex=1.5, col=txcol)
text(x=0.17+0.0065, y=0.63, labels= bquote(pi[1] == 0.1), cex=1.5, col=txcol)
text(x=0.17, y=0.56, labels= bquote(paste(pi,"'") == 0), cex=1.5, col=txcol)
text(x=0.17+0.0065, y=0.49, labels= bquote("p'" == 0.1), cex=1.5, col=txcol)
text(x=0.17+0.005, y=0.42, labels= bquote(gamma == 0.8), cex=1.5, col=txcol)


mtext(side=1, "Point Estimates", line=2)
mtext(side=2, "Cumulative Distiburion Function", line=2)
title("B", adj=0, cex.main=2, font=2, line=0.5)
box()

############################################################
# PANEL C
############################################################

pi=0.1
plot(test[[2]], type="n", ylim=c(0,0.4),
     xlab="", las=1,
     ylab="")
points(test[[2]], pch=16, col=scales::alpha("dimgray",1))
arrows(x0=1:length(test[[2]]), x1=1:length(test[[2]]),
       y0=test[[3]], y1=test[[4]],
       length=0, col=scales::alpha("gray60",0.9))
abline(h=0, col="dimgray", lwd=1.5)
abline(h=pi, col="firebrick4", lty=1, lwd=1.5)
legend(x=0,y=0.42, legend=c("Quantity of interest"),
       lty=c(1), lwd=c(1.5),
       col=c("firebrick4"),
       box.lty=0, cex=1.5)
text(x=18, y=0.33, labels= bquote("n" == 1000), cex=1.5, col=txcol)
text(x=15, y=0.30, labels= bquote(pi == 0.1), cex=1.5, col=txcol)
text(x=12.5, y=0.27, labels= bquote(paste(pi,"'") == 0), cex=1.5, col=txcol)
text(x=15, y=0.24, labels= bquote("p" == 0.1), cex=1.5, col=txcol)
text(x=16, y=0.21, labels= bquote("p'" == 0.1), cex=1.5, col=txcol)
text(x=14.5, y=0.18, labels= bquote(gamma == 0.8), cex=1.5, col=txcol)

mtext(side=1, "Simulation Index", line=2)
mtext(side=2, "Prevalence of Sensitive Attribute", line=2)
title("C", adj=0, cex.main=2, font=2, line=0.5)
box()

dev.off()


end <- Sys.time()
start - end

#############################################################################
# END OF THIS R SOURCE FILE
#############################################################################
