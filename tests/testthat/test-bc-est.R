test_that("bc_est returns the documented prevalence output", {
  fit <- bc_est(
    Y = Y, A = A, p = 0.15, p.prime = 0.15,
    data = cmdata, seed = 20260901L
  )

  expect_named(fit, c("Results", "Stats"))
  expect_equal(dim(fit$Results), c(2L, 4L))
  expect_equal(
    rownames(fit$Results),
    c("Naive Crosswise", "Bias-Corrected")
  )
  expect_true(all(fit$Results[, "Estimate"] >= 0 & fit$Results[, "Estimate"] <= 1))
  expect_equal(unname(fit$Stats[1, "Sample Size"]), nrow(cmdata))
})

test_that("bc_est treats unit weights like unweighted data", {
  unit_weight_data <- transform(cmdata, unit_weight = 1)

  unweighted <- bc_est(
    Y = Y, A = A, p = 0.15, p.prime = 0.15,
    data = unit_weight_data, seed = 20260901L
  )
  weighted <- bc_est(
    Y = Y, A = A, weight = unit_weight, p = 0.15, p.prime = 0.15,
    data = unit_weight_data, seed = 20260901L
  )

  expect_identical(weighted, unweighted)
})

test_that("bc_est removes incomplete survey rows before estimating", {
  incomplete_data <- cmdata
  incomplete_data$Y[1] <- NA_integer_
  incomplete_data$A[2] <- NA_integer_
  complete_data <- incomplete_data[complete.cases(incomplete_data), , drop = FALSE]

  from_incomplete_data <- bc_est(
    Y = Y, A = A, p = 0.15, p.prime = 0.15,
    data = incomplete_data, seed = 20260901L
  )
  from_complete_data <- bc_est(
    Y = Y, A = A, p = 0.15, p.prime = 0.15,
    data = complete_data, seed = 20260901L
  )

  expect_identical(from_incomplete_data, from_complete_data)
  expect_equal(
    unname(from_incomplete_data$Stats[1, "Sample Size"]),
    nrow(complete_data)
  )
})
