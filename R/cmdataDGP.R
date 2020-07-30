
i.logit <- function(XB){ exp(XB)/(1 + exp(XB))}
set.seed(1234)

N = 2000
beta  = c(-1.5, 0.5, 0.02)                           # Unknown parameter
theta = c(2, -0.1, -0.01)                            # Unknown parameter
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

sample.prob = NA
sample.prob = ifelse(female==1, 0.3, 0.7)
sweight = 1/sample.prob

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

dat <- cbind(Y, A, female, age, p, p2)                          # Observed data
cmdata2 <- as.data.frame(dat)
save(cmdata2, file="cmdata2.RData")



