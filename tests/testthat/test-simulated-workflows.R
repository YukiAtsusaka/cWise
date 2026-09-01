test_that("bc_est works with package example data", {
  set.seed(20260825)

  prevalence <- bc_est(Y = Y, A = A, p = 0.15, p.prime = 0.15, data = cmdata)
  expect_true(is.matrix(prevalence$Results))
  expect_true(all(is.finite(prevalence$Results)))
})

test_that("simulation and regression workflows work on simulated data", {
  set.seed(20260826)
  simulation <- sim_cwdata(
    N.sim = 2, sample = 80, prevalence = 0.15, p = 0.20, p.prime = 0.20,
    gamma = 0.80, direct = 0.08
  )
  expect_named(simulation, c("Results", "BiasCorrectEst", "BiasCorrectLow",
                             "BiasCorrectHigh", "EstimatedBias", "RelativeLengthCI"))

  graphics_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(graphics_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(sim_estimates(
    sample = 80, prevalence = 0.15, p = 0.20, p.prime = 0.20,
    gamma = 0.80, direct = 0.08, sim.results = simulation
  ))

  power <- sim_power(
    N.sim = 2, sample = 80, pi.null = 0.05, pi.alt = 0.15,
    p = 0.20, p.prime = 0.20, gamma = 0.80, direct = 0.08
  )
  expect_true(is.finite(power))

  sample_grid <- sim_power_N(
    N.sim = 1, prevalence = 0.15, p = 0.20, p.prime = 0.20,
    gamma = 0.80, direct = 0.08
  )
  expect_equal(nrow(sample_grid), 7L)

  outcome_fit <- cmreg(Y ~ female + age, anchor = A, p = 0.10, p.prime = 0.15,
                       data = cmdata2, n.start = 1L)
  outcome_prediction <- cmpredict(
    outcome_fit, newdata = data.frame(female = 1, age = 30), nsim = 20, seed = 1
  )
  expect_equal(nrow(outcome_prediction), 1L)

  predictor_fit <- cmreg_p(V ~ age + female, crosswise = Y, anchor = A,
                           p = 0.10, p.prime = 0.15,
                           data = cmdata3, n.start = 1L)
  predictor_prediction <- cmpredict_p(
    predictor_fit, newdata = data.frame(age = 30, female = 1), nsim = 20, seed = 1
  )
  expect_equal(nrow(predictor_prediction), 2L)
})
