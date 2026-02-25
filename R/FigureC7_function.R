#############################################################################
# FigureC7_function.R
# Aim: to offer two help functions to replicate Figure C7 in Atsusaka and Stevenson (2021)
# Run time: about 1 seconds
#############################################################################

###################################################################################
sim.curve <- function(N.sim, pi.null, pi.alt, p, p.prime, gamma, direct){

n_vec <- c(1, seq(from=500, to=2500, by=500))
out_n <- NA
out_power <- NA
alpha <- 0.05
mu0 <- pi.null
mu1 <- pi.alt

for(i in 1:length(n_vec)){

n = n_vec[i]        
print(paste0("Now, working on the sample size of ", n))

H0 <- sim.power(N.sim=N.sim, sample=n, pi=pi.null, p=p, p.prime= p.prime, gamma=gamma, direct=direct)
H1 <- sim.power(N.sim=N.sim, sample=n, pi=pi.alt,  p=p, p.prime= p.prime, gamma=gamma, direct=direct)

# POWER 
# P("H1"|H1) = Phi( [mu1 - mu0 + z_alpha * sigma0]/sigma1), where
# Phi: normal CDF
# mu0: mean under H0
# mu1: mean under H1
# sigma0: simulated sigma under H0
# sigma1: simulated sigma under H1

sim.sigma0 <- sd(H0$BiasCorrectEst) # Sample Standard Deviation from Simulated Data
sim.sigma1 <- sd(H1$BiasCorrectEst) # Sample Standard Deviation from Simulated Data
c_one_minus_alpha <- quantile(H0$BiasCorrectEst, 1-alpha)
z_one_minus_alpha <- (c_one_minus_alpha - mu0)/sim.sigma0
z_alpha <- -1*z_one_minus_alpha     # Flipping the sign

A <- (mu1 - mu0 + z_alpha*sim.sigma0)/sim.sigma1
power <- pnorm(A, mean=0, sd=1) # Probability of rejecting H0 given H1 is right = Power

out_power[i] <- power
out_n[i] <- n
} # END OF THE LOOP (i)


out <- list(simulated_power = out_power,
            sample_size = out_n)

} # END OF THE FUNCTION
###################################################################################

#############################################################################
# (1) FIND ROUGHT IDEA ABOUT SAMPLE SIZE
#############################################################################

sim.power.N <- function(N.sim=50, pi, p, p.prime, gamma, direct){

# RECODING AND ADDING FUNCTIONALITIES
roundy <- function(x){round(x, digits=3)}

p2=p.prime
Nvec <- c(100,500,1000,1500,2000,2500,3000)
pb <- txtProgressBar(min=1, max=length(Nvec), initial=0, style=3) # PROGRESS BAR  

coverage.zero.save <- NA
coverage.direct.save <- NA


for(n in 1:length(Nvec)){

setTxtProgressBar(pb,n) # SET PROGRESS BAR     
N = Nvec[n]   # SAMPLE SIZE 
  

# STORAGE OF SIMULATION RESULTS
naive.cover <- NA
bc.cover <- NA
naive.pred <- NA
bc.pred <- NA
bc.high.save <- NA
bc.low.save <- NA
REbc <- NA
est.bias <- NA
bc.cover0 <- NA
bc.cover.direct <- NA

#############################################################
# REPREAT THE SIMULATION N.sim TIMES WITH THE SAME PARAMETERS
#############################################################
for(i in 1:N.sim){  
  
  Attentive = rbinom(n=N, size=1, prob=gamma) # Attentive R or not
  
  # SENSITIVE QUESTION OF INTERST
  StatementA = rbinom(n=N, size=1, prob=pi)           # Sesitive item
  StatementB = rbinom(n=N, size=1, prob=p)            # Randomization item
  True.res = ifelse(StatementA != StatementB, 0, 1)   # Survey answer 
  
  Obs.res = rep(NA, N)
  Obs.res[Attentive==1] = True.res[Attentive==1]      # Observed answer for attentive Rs
  Obs.res[Attentive==0] = rbinom(n=length(N[Attentive==0]),
                                 size=1, prob=0.5)           # Observed answer, otherwise
  
  # NON-SENSITIVE AUXIALLY QUESTION
  StatementC = 0                                      # Attention "c"heck item (True prevalence is 0 for everyone)
  StatementD = rbinom(n=N, size=1, prob=p2)           # Anchor question
  True.res.prime = ifelse(StatementC != StatementD, 0, 1)   # Survey answer in Anchor Question
  
  Obs.res.prime = rep(NA, N)
  Obs.res.prime[Attentive==1] = True.res.prime[Attentive==1] # Observed answer for attentive Rs
  Obs.res.prime[Attentive==0] = rbinom(n=length(N[Attentive==0]),
                                       size=1, prob=0.5)            # Observed answer, otherwise
  
  dat <- cbind(Obs.res, True.res, Attentive, Obs.res.prime, True.res.prime)
  dat <- as_tibble(dat)
  head(dat) # See the D
  
  
  # (2) NAIVE CROSSWISE MODEL
  lambda.hat = mean(Obs.res) # Observed proportion of YESYES or NONO
  pi.hat.naive = (lambda.hat+p-1)/(2*p-1)
  
  pi.hat.naive.var = (pi.hat.naive*(1-pi.hat.naive))/(N-1) +
    (p*(1-p))/((N-1)*((2*p-1)^2))
  pi.hat.naive.sd = sqrt(pi.hat.naive.var)     
  
  naive.low  = pi.hat.naive-1.96*pi.hat.naive.sd; naive.low  = ifelse(naive.low > 0,  naive.low, 0)
  naive.high = pi.hat.naive+1.96*pi.hat.naive.sd; naive.high = ifelse(naive.high < 1, naive.high, 1)
  
  
  # (3) BIAS CORRECTED CROSSWISE MODEL
  gamma.hat = (mean(Obs.res.prime)-0.5)/(0.5-p2) # Estimated level of inattentiveness
  Bias.hat = (1/2)*((lambda.hat-0.5)/(p-0.5)) - (1/(2*gamma.hat))*((lambda.hat-0.5)/(p-0.5))
  Bias.hat
  
  pi.hat.bc = pi.hat.naive - Bias.hat # Bias Correction
  pi.hat.bc <- ifelse(pi.hat.bc < 0, 0, ifelse(pi.hat.bc > 1, 1, pi.hat.bc))
  
# BOOTSTRAPPING (200 TIMES)
  bs <- NA
  for(b in 1:500){
    index <- sample(1:nrow(dat), size=dim(dat), replace=TRUE)   # RESAMPLE WITH REPLACEMENT
    bs.dat <- dat[index,]                                       # BOOTSTRAPPED DATA
    
    bs.lambda.hat = mean(bs.dat$Obs.res)                    # Observed proportion of YESYES or NONO
    bs.pi.hat.naive = (bs.lambda.hat+p-1)/(2*p-1)
    bs.gamma.hat = (mean(bs.dat$Obs.res.prime)-0.5)/(0.5-p2) # Estimated level of inattentiveness
    bs.bias.hat = (1/2)*((bs.lambda.hat-0.5)/(p-0.5)) - (1/(2*bs.gamma.hat))*((bs.lambda.hat-0.5)/(p-0.5)) 
    
    bs[b] = bs.pi.hat.naive - bs.bias.hat # Bias Correction within Bootstrapping
    
    bs[b] <- ifelse(bs[b] < 0, 0, bs[b])
    bs[b] <- ifelse(bs[b] > 1, 1, bs[b])
  } # END OF BOOSTRAPPING

    
  pi.hat.bc.var = var(bs)
  pi.hat.bc.sd = sqrt(pi.hat.bc.var)
  
  bc.low  = quantile(bs, prob=0.025, na.rm=T)
  bc.low  = ifelse(bc.low > 0,  bc.low, 0)
  bc.high = quantile(bs, prob=0.975, na.rm=T)
  bc.high = ifelse(bc.high < 1, bc.high, 1)
  
  
# STORING RESULTS
  est.bias[i] <- Bias.hat
  
  naive.pred[i] <- pi.hat.naive
  bc.pred[i] <- pi.hat.bc
  bc.high.save[i] <- bc.high
  bc.low.save[i] <- bc.low  

  naive.cover[i] <- (naive.low <= pi & pi <= naive.high) # 95CI contains pi
  bc.cover[i] <-    (bc.low <= pi & pi <= bc.high)       # 95CI contains pi

  bc.cover0[i] <- (bc.low <= 0)                                 # 95%CI contains 0
  bc.cover.direct[i] <- (bc.low <= direct & direct <= bc.high ) # 95%CI contains Direct Est  

  REbc[i] <- (bc.high-bc.low)/(naive.high-naive.low) # LENGTH OF 95CI
# END OF ONE SIMULATION
}
#############################################################
# END OF ALL SIMULATIONS for one n
#############################################################

#############################################################
# ORGANIZING OUTPUTS
#############################################################

coverage.zero <- mean(bc.cover0, na.rm=T)*100
coverage.direct <- mean(bc.cover.direct, na.rm=T)*100

coverage.zero.save[n] <- coverage.zero     # SAVE FOR EACH n
coverage.direct.save[n] <- coverage.direct # SAVE FOR EACH n
} # END OF OUTER LOOP (n)


#########
# OUTPUT
#########

out <- data.frame(CoverageZero = coverage.zero.save,
                  CoverageDirect = coverage.direct.save)


return(out)  
} # END OF THIS FUNCTION


#############################################################################
# (2) VISUALIZE POSSIBLE ESTIMATES OF THE BIAS-CORRECTED ESTIMATOR
#############################################################################

sim.power <- function(N.sim=500, sample, pi, p, p.prime, gamma, direct){

# RECODING AND ADDING FUNCTIONALITIES
roundy <- function(x){round(x, digits=3)}
N=sample
p2=p.prime
pb <- txtProgressBar(min=1, max=N.sim, initial=0, style=3) # PROGRESS BAR  

# STORAGE OF SIMULATION RESULTS
naive.cover <- NA
bc.cover <- NA
naive.pred <- NA
bc.pred <- NA
bc.high.save <- NA
bc.low.save <- NA
REbc <- NA
est.bias <- NA
bc.cover0 <- NA
bc.cover.direct <- NA

#############################################################
# REPREAT THE SIMULATION N.sim TIMES WITH THE SAME PARAMETERS
#############################################################
for(i in 1:N.sim){  
  
  setTxtProgressBar(pb,i) # SET PROGRESS BAR 
  
  Attentive = rbinom(n=N, size=1, prob=gamma)# Attentive R or not
  
  # SENSITIVE QUESTION OF INTERST
  StatementA = rbinom(n=N, size=1, prob=pi)           # Sesitive item
  StatementB = rbinom(n=N, size=1, prob=p)            # Randomization item
  True.res = ifelse(StatementA != StatementB, 0, 1)   # Survey answer 
  
  Obs.res = rep(NA, N)
  Obs.res[Attentive==1] = True.res[Attentive==1]      # Observed answer for attentive Rs
  Obs.res[Attentive==0] = rbinom(n=length(N[Attentive==0]),
                                 size=1, prob=0.5)           # Observed answer, otherwise
  
  # NON-SENSITIVE AUXIALLY QUESTION
  StatementC = 0                                      # Attention "c"heck item (True prevalence is 0 for everyone)
  StatementD = rbinom(n=N, size=1, prob=p2)           # Anchor question
  True.res.prime = ifelse(StatementC != StatementD, 0, 1)   # Survey answer in Anchor Question
  
  Obs.res.prime = rep(NA, N)
  Obs.res.prime[Attentive==1] = True.res.prime[Attentive==1] # Observed answer for attentive Rs
  Obs.res.prime[Attentive==0] = rbinom(n=length(N[Attentive==0]),
                                       size=1, prob=0.5)            # Observed answer, otherwise
  
  dat <- cbind(Obs.res, True.res, Attentive, Obs.res.prime, True.res.prime)
  dat <- as_tibble(dat)
  head(dat) # See the D
  
  
  # (2) NAIVE CROSSWISE MODEL
  lambda.hat = mean(Obs.res) # Observed proportion of YESYES or NONO
  pi.hat.naive = (lambda.hat+p-1)/(2*p-1)
  
  pi.hat.naive.var = (pi.hat.naive*(1-pi.hat.naive))/(N-1) +
    (p*(1-p))/((N-1)*((2*p-1)^2))
  pi.hat.naive.sd = sqrt(pi.hat.naive.var)     
  
  naive.low  = pi.hat.naive-1.96*pi.hat.naive.sd; naive.low  = ifelse(naive.low > 0,  naive.low, 0)
  naive.high = pi.hat.naive+1.96*pi.hat.naive.sd; naive.high = ifelse(naive.high < 1, naive.high, 1)
  
  
  # (3) BIAS CORRECTED CROSSWISE MODEL
  gamma.hat = (mean(Obs.res.prime)-0.5)/(0.5-p2) # Estimated level of inattentiveness
  Bias.hat = (1/2)*((lambda.hat-0.5)/(p-0.5)) - (1/(2*gamma.hat))*((lambda.hat-0.5)/(p-0.5))
  Bias.hat
  
  pi.hat.bc = pi.hat.naive - Bias.hat # Bias Correction
  pi.hat.bc <- ifelse(pi.hat.bc < 0, 0, ifelse(pi.hat.bc > 1, 1, pi.hat.bc))
  
# BOOTSTRAPPING (200 TIMES)
  bs <- NA
  for(b in 1:500){
    index <- sample(1:nrow(dat), size=dim(dat), replace=TRUE)   # RESAMPLE WITH REPLACEMENT
    bs.dat <- dat[index,]                                       # BOOTSTRAPPED DATA
    
    bs.lambda.hat = mean(bs.dat$Obs.res)                    # Observed proportion of YESYES or NONO
    bs.pi.hat.naive = (bs.lambda.hat+p-1)/(2*p-1)
    bs.gamma.hat = (mean(bs.dat$Obs.res.prime)-0.5)/(0.5-p2) # Estimated level of inattentiveness
    bs.bias.hat = (1/2)*((bs.lambda.hat-0.5)/(p-0.5)) - (1/(2*bs.gamma.hat))*((bs.lambda.hat-0.5)/(p-0.5)) 
    
    bs[b] = bs.pi.hat.naive - bs.bias.hat # Bias Correction within Bootstrapping
    
    bs[b] <- ifelse(bs[b] < 0, 0, bs[b])
    bs[b] <- ifelse(bs[b] > 1, 1, bs[b])
  } # END OF BOOSTRAPPING

    
  pi.hat.bc.var = var(bs)
  pi.hat.bc.sd = sqrt(pi.hat.bc.var)
  
  bc.low  = quantile(bs, prob=0.025, na.rm=T)
  bc.low  = ifelse(bc.low > 0,  bc.low, 0)
  bc.high = quantile(bs, prob=0.975, na.rm=T)
  bc.high = ifelse(bc.high < 1, bc.high, 1)
  
  
# STORING RESULTS
  est.bias[i] <- Bias.hat
  
  naive.pred[i] <- pi.hat.naive
  bc.pred[i] <- pi.hat.bc
  bc.high.save[i] <- bc.high
  bc.low.save[i] <- bc.low  

  naive.cover[i] <- (naive.low <= pi & pi <= naive.high) # 95CI contains pi
  bc.cover[i] <-    (bc.low <= pi & pi <= bc.high)       # 95CI contains pi

  bc.cover0[i] <- (bc.low <= 0)                                 # 95%CI contains 0
  bc.cover.direct[i] <- (bc.low <= direct & direct <= bc.high ) # 95%CI contains Direct Est  

  REbc[i] <- (bc.high-bc.low)/(naive.high-naive.low) # LENGTH OF 95CI
# END OF ONE SIMULATION
}
#############################################################
# END OF ALL SIMULATIONS
#############################################################

#############################################################
# ORGANIZING OUTPUTS
#############################################################
  
# Bias: sum_{i=1}^{n}(\hat\pi - \pi)/n
bias.naive.cm      <- sum(naive.pred - pi)/N.sim 
bias.bias.correct  <- sum(bc.pred - pi)/N.sim

# RMSE: \sqrt( (sum_{i=1}^{n}(\hat\pi - \pi)^2)/n )
rmse.naive.cm      <- sum( (naive.pred - pi)^2 )/N.sim
rmse.bias.correct  <- sum( (bc.pred - pi)^2 )/N.sim

# What Percentage of Times Did 95CIs Capture the Ground Truth?
coverage.naive.cm     <- mean(naive.cover, na.rm =T)*100
coverage.bias.correct <- mean(bc.cover, na.rm=T)*100

coverage.zero <- mean(bc.cover0, na.rm=T)*100
coverage.direct <- mean(bc.cover.direct, na.rm=T)*100

# SORTING
bc.pred.sort <- sort(bc.pred)
bc.low.save.sort <- bc.low.save[order(bc.pred)]  
bc.high.save.sort <- bc.high.save[order(bc.pred)]  

#########
# OUTPUT
#########
out <- list()

outvec <- c(bias.naive.cm, bias.bias.correct,
            rmse.naive.cm, rmse.bias.correct,
            coverage.naive.cm, coverage.bias.correct,
            coverage.zero, coverage.direct)
outvec <- sapply(outvec, roundy)
names(outvec) <- c("Average Bias (naive)", "Average Bias (bias-corrected)",
                        "RMSE (naive)", "RMSE (bias-corrected)",
                        "% Coverage (naive)", "% Coverage (bias-corrected)",
                   "% Coverage (bias-corrected) of 0",
                   "% Coverage (bias-corrected) of Direct Estimate")

out <- list(Results = outvec,
            BiasCorrectEst = bc.pred.sort,       # POINT ESTIMATES FROM THE BIAS-CORRECTED ESTIMATOR
            BiasCorrectLow = bc.low.save.sort,   # LOWER BOUND OF THE 95%CI
            BiasCorrectHigh = bc.high.save.sort, # UPPER BOUND OF THE 95%CI
            EstimatedBias = est.bias, # ESTIMATD BIAS
            RelativeLengthCI = REbc)  # RELATIVE LENGTH OF 95% CI


return(out)  
}



#############################################################################
# END OF THIS R SOURCE FILE
#############################################################################