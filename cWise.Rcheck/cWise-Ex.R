pkgname <- "cWise"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "cWise-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('cWise')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bc_est")
### * bc_est

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bc_est
### Title: bc_est
### Aliases: bc_est bc.est

### ** Examples

bc_est(Y=Y, A=A, p=0.15, p.prime=0.15, data=cmdata)
bc_est(Y=Y, A=A, weight=weight, p=0.15, p.prime=0.15, data=cmdata)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bc_est", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cmBound")
### * cmBound

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cmBound
### Title: cmBound
### Aliases: cmBound

### ** Examples

p <- cmBound(lambda.hat=0.6385, p=0.25, N=310, dq=0.073, N.dq=310)
p

p <- p + ggplot2::ggtitle("Sensitivity Analysis") +
         ggplot2::theme(plot.title = ggplot2::element_text(size=20, face="bold"))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cmBound", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cmdata")
### * cmdata

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cmdata
### Title: Simulated Crosswise Survey Data I
### Aliases: cmdata
### Keywords: datasets

### ** Examples


data(cmdata)
head(cmdata)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cmdata", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cmdata2")
### * cmdata2

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cmdata2
### Title: Simulated Crosswise Survey Data II
### Aliases: cmdata2
### Keywords: datasets

### ** Examples


data(cmdata2)
head(cmdata2)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cmdata2", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cmdata3")
### * cmdata3

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cmdata3
### Title: Simulated Crosswise Survey Data III
### Aliases: cmdata3
### Keywords: datasets

### ** Examples


data(cmdata3)
head(cmdata3)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cmdata3", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cmpredict")
### * cmpredict

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cmpredict
### Title: cmpredict
### Aliases: cmpredict

### ** Examples

m <- cmreg(Y ~ female + age, anchor = A, p = 0.1, p.prime = 0.15,
           data = cmdata2)
predictions <- cmpredict(m, typical = c(age = 30),
                          zval = c(female = 0, female = 1))
predictions



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cmpredict", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cmpredict_p")
### * cmpredict_p

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cmpredict_p
### Title: cmpredict_p
### Aliases: cmpredict_p cmpredict.p

### ** Examples

m2 <- cmreg_p(V ~ age + female, crosswise = Y, anchor = A, p = 0.1,
              p.prime = 0.15, data = cmdata3)
predictions <- cmpredict_p(m2, newdata = data.frame(age = 30, female = 1))
predictions



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cmpredict_p", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cmreg")
### * cmreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cmreg
### Title: cmreg
### Aliases: cmreg

### ** Examples

m <- cmreg(Y ~ female + age, anchor = A, p = 0.1, p.prime = 0.15,
           data = cmdata2)
m



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cmreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cmreg_p")
### * cmreg_p

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cmreg_p
### Title: cmreg_p
### Aliases: cmreg_p cmreg.p

### ** Examples

m2 <- cmreg_p(V ~ age + female, crosswise = Y, anchor = A, p = 0.1,
              p.prime = 0.15, data = cmdata3)
m2



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cmreg_p", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim_cwdata")
### * sim_cwdata

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim_cwdata
### Title: Simulate Survey Data and Compute Bias-Corrected Estimates
### Aliases: sim_cwdata sim.cwdata

### ** Examples

## Not run: 
##D # Basic usage
##D result <- sim_cwdata(
##D   N.sim = 100,
##D   sample = 500,
##D   prevalence = 0.1,
##D   p = 0.1,
##D   p.prime = 0.1,
##D   gamma = 0.8,
##D   direct = 0.05
##D )
##D print(result$Results)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim_cwdata", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim_estimates")
### * sim_estimates

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim_estimates
### Title: Simulate and Plot Panel C of Figure C7
### Aliases: sim_estimates sim.estimates

### ** Examples

## Not run: 
##D # Replicate Panel C of Figure C7
##D sim_estimates(
##D   N.sim   = 100,
##D   sample  = 1000,
##D   prevalence = 0.1,
##D   p       = 0.1,
##D   p.prime = 0.1,
##D   gamma   = 0.8,
##D   direct  = 0.1
##D )
##D 
##D # Re-use pre-computed simulation results
##D res <- sim_cwdata(N.sim = 100, sample = 1000, prevalence = 0.1,
##D                   p = 0.1, p.prime = 0.1, gamma = 0.8, direct = 0.1)
##D sim_estimates(sample = 1000, prevalence = 0.1, p = 0.1, p.prime = 0.1,
##D               gamma = 0.8, direct = 0.1, sim.results = res)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim_estimates", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim_power")
### * sim_power

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim_power
### Title: Compute Statistical Power for a Fixed Sample Size
### Aliases: sim_power sim.power

### ** Examples

# Compute power at a fixed sample size of 1000
## Not run: 
##D pwr <- sim_power(
##D   N.sim  = 500,
##D   sample = 1000,
##D   pi.null = 0,
##D   pi.alt  = 0.1,
##D   p       = 0.1,
##D   p.prime = 0.1,
##D   gamma   = 0.8,
##D   direct  = 0.02
##D )
##D cat(sprintf("Estimated power: %.3f\n", pwr))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim_power", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sim_power_N")
### * sim_power_N

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sim_power_N
### Title: Determine Sample Size for Desired Coverage Properties
### Aliases: sim_power_N sim.power.N

### ** Examples

# Find sample size needed to reliably exclude zero
## Not run: 
##D result <- sim_power_N(
##D   N.sim = 50,
##D   prevalence = 0.1,
##D   p = 0.1,
##D   p.prime = 0.1,
##D   gamma = 0.8,
##D   direct = 0.05
##D )
##D 
##D print(result)
##D 
##D # Visualize results
##D plot(result$SampleSize, result$CoverageZero,
##D      type = "b", xlab = "Sample Size",
##D      ylab = "% of CIs Including Zero",
##D      main = "Precision vs Sample Size")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sim_power_N", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
