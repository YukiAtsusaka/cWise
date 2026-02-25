#############################################################################
# power_analysis.R
# Aim: Power analysis functions for crosswise model with bias correction
# Based on: Appendix C5 methodology
# Fixed version of FigureC7_function.R with all errors corrected
#############################################################################

###################################################################################
# sim.curve: Generate power curves for different sample sizes
#
# This function computes statistical power across a range of sample sizes
# for a one-sided hypothesis test of H0: π ≤ π0 vs H1: π > π0
#
# Parameters:
#   N.sim: Number of Monte Carlo simulations per sample size
#   pi.null: Prevalence rate under the null hypothesis (π0)
#   pi.alt: True prevalence rate under the alternative hypothesis (π1)
#   p: Probability for the randomization item in the sensitive question
#   p.prime: Probability for the anchor question (non-sensitive)
#   gamma: Proportion of attentive respondents
#   direct: Direct questioning estimate (for coverage comparison)
#
# Returns:
#   List with simulated_power and sample_size vectors
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

    H0 <- sim.power(N.sim=N.sim, sample=n, pi=pi.null, p=p, p.prime=p.prime, gamma=gamma, direct=direct)
    H1 <- sim.power(N.sim=N.sim, sample=n, pi=pi.alt,  p=p, p.prime=p.prime, gamma=gamma, direct=direct)

    # POWER CALCULATION
    # Power = P(Reject H0 | H1 is true) = Φ((μ1 - μ0 - c·σ0) / σ1)
    # where Φ is the standard normal CDF

    sim.sigma0 <- sd(H0$BiasCorrectEst) # Sample SD from simulated data under H0
    sim.sigma1 <- sd(H1$BiasCorrectEst) # Sample SD from simulated data under H1
    c_one_minus_alpha <- quantile(H0$BiasCorrectEst, 1-alpha)
    z_one_minus_alpha <- (c_one_minus_alpha - mu0)/sim.sigma0
    z_alpha <- -1*z_one_minus_alpha     # Flipping the sign

    A <- (mu1 - mu0 + z_alpha*sim.sigma0)/sim.sigma1
    power <- pnorm(A, mean=0, sd=1) # Probability of rejecting H0 given H1 is true

    out_power[i] <- power
    out_n[i] <- n
  } # END OF LOOP (i)


  out <- list(simulated_power = out_power,
              sample_size = out_n)

  return(out)
} # END OF sim.curve


###################################################################################
# sim.power.N: Determine rough sample size needed for desired coverage
#
# This function tests multiple sample sizes to see which achieves desired
# coverage properties (e.g., 95% CI excludes 0 or includes direct estimate)
#
# Parameters:
#   N.sim: Number of Monte Carlo simulations (default: 50)
#   pi: True prevalence rate of the sensitive attribute
#   p: Probability for the randomization item
#   p.prime: Probability for the anchor question
#   gamma: Proportion of attentive respondents
#   direct: Direct questioning estimate
#
# Returns:
#   Data frame with coverage statistics for different sample sizes
###################################################################################
sim.power.N <- function(N.sim=50, pi, p, p.prime, gamma, direct){

  # Helper function for rounding
  roundy <- function(x){round(x, digits=3)}

  p2 <- p.prime
  Nvec <- c(100, 500, 1000, 1500, 2000, 2500, 3000)
  pb <- txtProgressBar(min=1, max=length(Nvec), initial=0, style=3)

  coverage.zero.save <- NA
  coverage.direct.save <- NA

  for(n in 1:length(Nvec)){

    setTxtProgressBar(pb, n)
    N <- Nvec[n]   # Sample size

    # Storage for simulation results
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
    # Repeat simulation N.sim times with the same parameters
    #############################################################
    for(i in 1:N.sim){

      Attentive <- rbinom(n=N, size=1, prob=gamma) # Attentive respondents

      # SENSITIVE QUESTION OF INTEREST
      StatementA <- rbinom(n=N, size=1, prob=pi)           # Sensitive item
      StatementB <- rbinom(n=N, size=1, prob=p)            # Randomization item
      True.res <- ifelse(StatementA != StatementB, 0, 1)   # True survey answer

      Obs.res <- rep(NA, N)
      Obs.res[Attentive==1] <- True.res[Attentive==1]      # Observed answer for attentive Rs
      # FIXED: Use sum(Attentive==0) instead of length(N[Attentive==0])
      Obs.res[Attentive==0] <- rbinom(n=sum(Attentive==0),
                                     size=1, prob=0.5)      # Random response for inattentive Rs

      # NON-SENSITIVE AUXILIARY QUESTION
      StatementC <- 0                                       # Attention check item (True prevalence is 0)
      StatementD <- rbinom(n=N, size=1, prob=p2)            # Anchor question
      True.res.prime <- ifelse(StatementC != StatementD, 0, 1)

      Obs.res.prime <- rep(NA, N)
      Obs.res.prime[Attentive==1] <- True.res.prime[Attentive==1]
      # FIXED: Use sum(Attentive==0) instead of length(N[Attentive==0])
      Obs.res.prime[Attentive==0] <- rbinom(n=sum(Attentive==0),
                                           size=1, prob=0.5)

      dat <- cbind(Obs.res, True.res, Attentive, Obs.res.prime, True.res.prime)
      dat <- as.data.frame(dat)  # FIXED: Use as.data.frame instead of as_tibble


      # (2) NAIVE CROSSWISE MODEL
      lambda.hat <- mean(Obs.res) # Observed proportion of same answers
      pi.hat.naive <- (lambda.hat + p - 1)/(2*p - 1)

      pi.hat.naive.var <- (pi.hat.naive*(1-pi.hat.naive))/(N-1) +
        (p*(1-p))/((N-1)*((2*p-1)^2))
      pi.hat.naive.sd <- sqrt(pi.hat.naive.var)

      naive.low  <- pi.hat.naive - 1.96*pi.hat.naive.sd
      naive.low  <- ifelse(naive.low > 0, naive.low, 0)
      naive.high <- pi.hat.naive + 1.96*pi.hat.naive.sd
      naive.high <- ifelse(naive.high < 1, naive.high, 1)


      # (3) BIAS CORRECTED CROSSWISE MODEL
      gamma.hat <- (mean(Obs.res.prime) - 0.5)/(0.5 - p2) # Estimated attention rate
      Bias.hat <- (1/2)*((lambda.hat - 0.5)/(p - 0.5)) - (1/(2*gamma.hat))*((lambda.hat - 0.5)/(p - 0.5))

      pi.hat.bc <- pi.hat.naive - Bias.hat # Bias correction
      pi.hat.bc <- ifelse(pi.hat.bc < 0, 0, ifelse(pi.hat.bc > 1, 1, pi.hat.bc))

      # BOOTSTRAPPING (500 iterations)
      bs <- NA
      for(b in 1:500){
        # FIXED: Use nrow(dat) instead of dim(dat)
        index <- sample(1:nrow(dat), size=nrow(dat), replace=TRUE)
        bs.dat <- dat[index, ]

        bs.lambda.hat <- mean(bs.dat$Obs.res)
        bs.pi.hat.naive <- (bs.lambda.hat + p - 1)/(2*p - 1)
        bs.gamma.hat <- (mean(bs.dat$Obs.res.prime) - 0.5)/(0.5 - p2)
        bs.bias.hat <- (1/2)*((bs.lambda.hat - 0.5)/(p - 0.5)) - (1/(2*bs.gamma.hat))*((bs.lambda.hat - 0.5)/(p - 0.5))

        bs[b] <- bs.pi.hat.naive - bs.bias.hat # Bias correction in bootstrap

        bs[b] <- ifelse(bs[b] < 0, 0, bs[b])
        bs[b] <- ifelse(bs[b] > 1, 1, bs[b])
      } # END OF BOOTSTRAPPING

      pi.hat.bc.var <- var(bs)
      pi.hat.bc.sd <- sqrt(pi.hat.bc.var)

      bc.low  <- quantile(bs, prob=0.025, na.rm=TRUE)
      bc.low  <- ifelse(bc.low > 0, bc.low, 0)
      bc.high <- quantile(bs, prob=0.975, na.rm=TRUE)
      bc.high <- ifelse(bc.high < 1, bc.high, 1)


      # STORING RESULTS
      est.bias[i] <- Bias.hat

      naive.pred[i] <- pi.hat.naive
      bc.pred[i] <- pi.hat.bc
      bc.high.save[i] <- bc.high
      bc.low.save[i] <- bc.low

      naive.cover[i] <- (naive.low <= pi & pi <= naive.high)
      bc.cover[i] <- (bc.low <= pi & pi <= bc.high)

      bc.cover0[i] <- (bc.low <= 0)
      bc.cover.direct[i] <- (bc.low <= direct & direct <= bc.high)

      REbc[i] <- (bc.high - bc.low)/(naive.high - naive.low)
    } # END OF ONE SIMULATION

    #############################################################
    # ORGANIZING OUTPUTS
    #############################################################

    coverage.zero <- mean(bc.cover0, na.rm=TRUE)*100
    coverage.direct <- mean(bc.cover.direct, na.rm=TRUE)*100

    coverage.zero.save[n] <- coverage.zero
    coverage.direct.save[n] <- coverage.direct
  } # END OF OUTER LOOP (n)


  #########
  # OUTPUT
  #########

  out <- data.frame(SampleSize = Nvec,
                    CoverageZero = coverage.zero.save,
                    CoverageDirect = coverage.direct.save)

  return(out)
} # END OF sim.power.N


###################################################################################
# sim.power: Simulate survey data and compute bias-corrected estimates
#
# Core simulation function that generates crosswise model data with inattentive
# respondents and computes bias-corrected prevalence estimates with bootstrap CIs
#
# Parameters:
#   N.sim: Number of Monte Carlo simulations
#   sample: Sample size per simulation
#   pi: True prevalence rate of the sensitive attribute
#   p: Probability for the randomization item
#   p.prime: Probability for the anchor question
#   gamma: Proportion of attentive respondents
#   direct: Direct questioning estimate (for comparison)
#
# Returns:
#   List containing:
#     - Results: Summary statistics (bias, RMSE, coverage)
#     - BiasCorrectEst: Vector of bias-corrected point estimates
#     - BiasCorrectLow/High: Bootstrap CI bounds
#     - EstimatedBias: Vector of estimated bias values
#     - RelativeLengthCI: Relative CI length vs naive estimator
###################################################################################
sim.power <- function(N.sim=500, sample, pi, p, p.prime, gamma, direct){

  # Helper function for rounding
  roundy <- function(x){round(x, digits=3)}

  N <- sample
  p2 <- p.prime
  pb <- txtProgressBar(min=1, max=N.sim, initial=0, style=3)

  # Storage for simulation results
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
  # Repeat simulation N.sim times with the same parameters
  #############################################################
  for(i in 1:N.sim){

    setTxtProgressBar(pb, i)

    Attentive <- rbinom(n=N, size=1, prob=gamma) # Attentive respondents

    # SENSITIVE QUESTION OF INTEREST
    StatementA <- rbinom(n=N, size=1, prob=pi)           # Sensitive item
    StatementB <- rbinom(n=N, size=1, prob=p)            # Randomization item
    True.res <- ifelse(StatementA != StatementB, 0, 1)   # True survey answer

    Obs.res <- rep(NA, N)
    Obs.res[Attentive==1] <- True.res[Attentive==1]      # Observed answer for attentive Rs
    # FIXED: Use sum(Attentive==0) instead of length(N[Attentive==0])
    Obs.res[Attentive==0] <- rbinom(n=sum(Attentive==0),
                                   size=1, prob=0.5)      # Random response for inattentive Rs

    # NON-SENSITIVE AUXILIARY QUESTION
    StatementC <- 0                                       # Attention check item
    StatementD <- rbinom(n=N, size=1, prob=p2)            # Anchor question
    True.res.prime <- ifelse(StatementC != StatementD, 0, 1)

    Obs.res.prime <- rep(NA, N)
    Obs.res.prime[Attentive==1] <- True.res.prime[Attentive==1]
    # FIXED: Use sum(Attentive==0) instead of length(N[Attentive==0])
    Obs.res.prime[Attentive==0] <- rbinom(n=sum(Attentive==0),
                                         size=1, prob=0.5)

    dat <- cbind(Obs.res, True.res, Attentive, Obs.res.prime, True.res.prime)
    dat <- as.data.frame(dat)  # FIXED: Use as.data.frame instead of as_tibble


    # (2) NAIVE CROSSWISE MODEL
    lambda.hat <- mean(Obs.res) # Observed proportion
    pi.hat.naive <- (lambda.hat + p - 1)/(2*p - 1)

    pi.hat.naive.var <- (pi.hat.naive*(1-pi.hat.naive))/(N-1) +
      (p*(1-p))/((N-1)*((2*p-1)^2))
    pi.hat.naive.sd <- sqrt(pi.hat.naive.var)

    naive.low  <- pi.hat.naive - 1.96*pi.hat.naive.sd
    naive.low  <- ifelse(naive.low > 0, naive.low, 0)
    naive.high <- pi.hat.naive + 1.96*pi.hat.naive.sd
    naive.high <- ifelse(naive.high < 1, naive.high, 1)


    # (3) BIAS CORRECTED CROSSWISE MODEL
    gamma.hat <- (mean(Obs.res.prime) - 0.5)/(0.5 - p2) # Estimated attention rate
    Bias.hat <- (1/2)*((lambda.hat - 0.5)/(p - 0.5)) - (1/(2*gamma.hat))*((lambda.hat - 0.5)/(p - 0.5))

    pi.hat.bc <- pi.hat.naive - Bias.hat # Bias correction
    pi.hat.bc <- ifelse(pi.hat.bc < 0, 0, ifelse(pi.hat.bc > 1, 1, pi.hat.bc))

    # BOOTSTRAPPING (500 iterations)
    bs <- NA
    for(b in 1:500){
      # FIXED: Use nrow(dat) instead of dim(dat)
      index <- sample(1:nrow(dat), size=nrow(dat), replace=TRUE)
      bs.dat <- dat[index, ]

      bs.lambda.hat <- mean(bs.dat$Obs.res)
      bs.pi.hat.naive <- (bs.lambda.hat + p - 1)/(2*p - 1)
      bs.gamma.hat <- (mean(bs.dat$Obs.res.prime) - 0.5)/(0.5 - p2)
      bs.bias.hat <- (1/2)*((bs.lambda.hat - 0.5)/(p - 0.5)) - (1/(2*bs.gamma.hat))*((bs.lambda.hat - 0.5)/(p - 0.5))

      bs[b] <- bs.pi.hat.naive - bs.bias.hat

      bs[b] <- ifelse(bs[b] < 0, 0, bs[b])
      bs[b] <- ifelse(bs[b] > 1, 1, bs[b])
    } # END OF BOOTSTRAPPING

    pi.hat.bc.var <- var(bs)
    pi.hat.bc.sd <- sqrt(pi.hat.bc.var)

    bc.low  <- quantile(bs, prob=0.025, na.rm=TRUE)
    bc.low  <- ifelse(bc.low > 0, bc.low, 0)
    bc.high <- quantile(bs, prob=0.975, na.rm=TRUE)
    bc.high <- ifelse(bc.high < 1, bc.high, 1)


    # STORING RESULTS
    est.bias[i] <- Bias.hat

    naive.pred[i] <- pi.hat.naive
    bc.pred[i] <- pi.hat.bc
    bc.high.save[i] <- bc.high
    bc.low.save[i] <- bc.low

    naive.cover[i] <- (naive.low <= pi & pi <= naive.high)
    bc.cover[i] <- (bc.low <= pi & pi <= bc.high)

    bc.cover0[i] <- (bc.low <= 0)
    bc.cover.direct[i] <- (bc.low <= direct & direct <= bc.high)

    REbc[i] <- (bc.high - bc.low)/(naive.high - naive.low)
  } # END OF ONE SIMULATION
  #############################################################
  # END OF ALL SIMULATIONS
  #############################################################

  #############################################################
  # ORGANIZING OUTPUTS
  #############################################################

  # Bias: sum(estimate - true_value)/N.sim
  bias.naive.cm <- sum(naive.pred - pi)/N.sim
  bias.bias.correct <- sum(bc.pred - pi)/N.sim

  # RMSE: sum((estimate - true_value)^2)/N.sim
  rmse.naive.cm <- sum((naive.pred - pi)^2)/N.sim
  rmse.bias.correct <- sum((bc.pred - pi)^2)/N.sim

  # Coverage: Percentage of times 95% CI captures the true value
  coverage.naive.cm <- mean(naive.cover, na.rm=TRUE)*100
  coverage.bias.correct <- mean(bc.cover, na.rm=TRUE)*100

  coverage.zero <- mean(bc.cover0, na.rm=TRUE)*100
  coverage.direct <- mean(bc.cover.direct, na.rm=TRUE)*100

  # Sort estimates for visualization
  bc.pred.sort <- sort(bc.pred)
  bc.low.save.sort <- bc.low.save[order(bc.pred)]
  bc.high.save.sort <- bc.high.save[order(bc.pred)]

  #########
  # OUTPUT
  #########
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
              BiasCorrectEst = bc.pred.sort,       # Point estimates
              BiasCorrectLow = bc.low.save.sort,   # Lower bound of 95% CI
              BiasCorrectHigh = bc.high.save.sort, # Upper bound of 95% CI
              EstimatedBias = est.bias,            # Estimated bias
              RelativeLengthCI = REbc)             # Relative length of 95% CI

  return(out)
}


#############################################################################
# END OF THIS R SOURCE FILE
#############################################################################
