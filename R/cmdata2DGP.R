
i.logit <- function(XB){ exp(XB)/(1 + exp(XB))}


set.seed(1234)
N = 2000

beta  = c(-1.5, 0.5, 0.02)                           # Unknown parameter (associating covariates and sensitive trais)
theta = c(2, -0.1, -0.01)                            # Unknown parameter (associating covariates and attentiveness)
female = rbinom(n=N, size=1, prob=0.5)               # 50% Male, 50% Female
age = rpois(n=N, lambda=30)                          # Age with mean 30
x = as.matrix(data.frame(rep(1, N), female, age))

XB =  x %*% beta                                     # Linear aggregator with betas
XT =  x %*% theta                                    # Linear aggregator with thetas
pi.x    = i.logit(XB)                                # Pi|Male
gamma.x = i.logit(XT)                                # Gamma|Male

p = 0.1                                              # Randomization probability
p2 = 0.15                                            # Randomization probability for the anchor question
Attentive = rbinom(n=N, size=1, prob=gamma.x)        # Attentive R or not

# glm(pi.x~x, family=binomial)                       # CHECK
# glm(gamma.x~x, family=binomial)                    # CHECK


# SENSITIVE QUESTION OF INTERST
StatementA = rbinom(n=N, size=1, prob=pi.x)          # Sesitive item
StatementB = rbinom(n=N, size=1, prob=p)             # Randomization item
True.res = ifelse(StatementA != StatementB, 0, 1)    # Survey answer

Obs.res = rep(NA, N)
Obs.res[Attentive==1] = True.res[Attentive==1]       # Observed answer for attentive Rs
Obs.res[Attentive==0] = rbinom(n=length(N[Attentive==0]),
                               size=1, prob=0.5)     # Observed answer, otherwise


# NON-SENSITIVE ANCHOR QUESTION
StatementC = 0                                       # Attention "c"heck item (True prevalence is 0 for everyone)
StatementD = rbinom(n=N, size=1, prob=p2)            # Anchor question
True.res.prime = ifelse(StatementC != StatementD, 0, 1) # Survey answer in Anchor Question

Obs.res.prime = rep(NA, N)
Obs.res.prime[Attentive==1] = True.res.prime[Attentive==1] # Observed answer for attentive Rs
Obs.res.prime[Attentive==0] = rbinom(n=length(N[Attentive==0]),
                                     size=1, prob=0.5)# Observed answer, otherwise

Y = Obs.res
A = Obs.res.prime


# PARAMETERS associating covariates + sensitive latent trait with the OUTCOME variable (NORMAL)
gamma0 = 0
gamma1 = 0.3
gamma2 = 0.01
gamma.cm = 1
SD = 1

Z = StatementA                                    # Having a sensitive attribute

V = gamma0*rep(1,N) + gamma1*female + gamma2*age + gamma.cm*Z + rnorm(mean=0, sd=SD, n=N)   # Outcome
gamma <- c(gamma0, gamma1, gamma2, gamma.cm)


dat <- cbind(V, Y, female, age, A, p, p2)                          # Observed data
dat <- as.data.frame(dat) # Data Frame
head(dat)


cmdata3 <- as.data.frame(dat)
save(cmdata3, file="cmdata3.RData")



