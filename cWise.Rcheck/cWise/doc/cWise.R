## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(cWise)

## ----data---------------------------------------------------------------------
head(cmdata)

## ----prevalence---------------------------------------------------------------
estimate <- bc_est(Y = Y, A = A, p = 0.15, p.prime = 0.15, data = cmdata)
estimate

## ----bounds, fig.width = 7, fig.height = 4------------------------------------
cmBound(lambda.hat = mean(cmdata$Y), p = 0.15, N = nrow(cmdata))

## ----power, eval = FALSE------------------------------------------------------
#  simulation <- sim_cwdata(
#    N.sim = 200, sample = 500, prevalence = 0.10, p = 0.15, p.prime = 0.15,
#    gamma = 0.80, direct = 0.05
#  )
#  simulation$Results
#  
#  sim_power(
#    N.sim = 500, sample = 1000, pi.null = 0.05, pi.alt = 0.10,
#    p = 0.15, p.prime = 0.15, gamma = 0.80, direct = 0.05
#  )

