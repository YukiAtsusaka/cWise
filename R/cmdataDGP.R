
set.seed(12)
N = 2000        # Sample size
pi = 0.1        # True prevelance
p = 0.15        # Randomization probability
p2 = 0.15       # Anchor randomization probability
gamma = 0.8    # Attentive rate
Attentive = rbinom(n=N, size=1, prob=gamma)         # Attentive R or not

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


female = rbinom(n=N, size=1, prob=0.8) # 50% Male, 50% Female --> But oversample female
w = ifelse(female==1, 1/(0.8), 1/(0.2))

#For package
cmdata <- dat %>% as_tibble %>%
       mutate(Y = Obs.res, A = Obs.res.prime, weight=w,
               p = p, p.prime = p2) %>%
       select(Y, A, weight, p, p.prime) %>% as.data.frame()
save(cmdata, file="cmdata.RData")

