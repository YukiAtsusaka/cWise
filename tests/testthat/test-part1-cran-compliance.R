test_that("bc_est seeding is reproducible without changing the caller RNG", {
  set.seed(20260830)
  before <- .Random.seed

  first <- bc_est(Y = Y, A = A, p = 0.15, p.prime = 0.15,
                  data = cmdata, seed = 111L)
  after <- .Random.seed
  second <- bc_est(Y = Y, A = A, p = 0.15, p.prime = 0.15,
                   data = cmdata, seed = 111L)

  expect_identical(after, before)
  expect_equal(first, second)
})

test_that("simulation APIs use prevalence and can be quiet", {
  expect_false("pi" %in% names(formals(sim_cwdata)))
  expect_true("prevalence" %in% names(formals(sim_cwdata)))
  expect_true("verbose" %in% names(formals(sim_cwdata)))

  set.seed(20260831)
  quiet <- capture.output(
    result <- sim_cwdata(N.sim = 1L, sample = 25L, prevalence = 0.1,
                         p = 0.15, p.prime = 0.15, gamma = 0.8,
                         direct = 0.05, verbose = FALSE)
  )
  expect_length(quiet, 0L)
  expect_named(result, c(
    "Results", "BiasCorrectEst", "BiasCorrectLow", "BiasCorrectHigh",
    "EstimatedBias", "RelativeLengthCI"
  ))

  set.seed(20260831)
  expect_warning(
    legacy <- sim.cwdata(N.sim = 1L, sample = 25L, pi = 0.1,
                          p = 0.15, p.prime = 0.15, gamma = 0.8,
                          direct = 0.05, verbose = FALSE),
    "deprecated"
  )
  expect_equal(legacy$Results, result$Results)
})

test_that("cmBound examples use current ggplot2 syntax", {
  expect_s3_class(
    cmBound(lambda.hat = 0.6385, p = 0.25, N = 310,
            dq = 0.073, N.dq = 310),
    "ggplot"
  )
})
