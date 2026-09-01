test_that("cmBound produces sensitivity plots", {
  expect_s3_class(
    cmBound(lambda.hat = 0.6385, p = 0.25, N = 310),
    "ggplot"
  )
  expect_s3_class(
    cmBound(lambda.hat = 0.6385, p = 0.25, N = 310,
            dq = 0.073, N.dq = 310),
    "ggplot"
  )
})

test_that("cmBound rejects undefined, invalid, or obsolete inputs", {
  expect_error(
    cmBound(lambda.hat = 0.6385, p = 0.25, N = 310, kappa = 0.5),
    "unused argument"
  )
  expect_error(
    cmBound(lambda.hat = 0.6385, p = 0.5, N = 310),
    "p.*0\\.5"
  )
  expect_error(
    cmBound(lambda.hat = 0.6385, p = 0, N = 310),
    "single finite probability"
  )
  expect_error(
    cmBound(lambda.hat = 0.6385, p = 1, N = 310),
    "single finite probability"
  )
  expect_error(
    cmBound(lambda.hat = 0.6385, p = 0.25, N = 310,
            dq = -0.01, N.dq = 310),
    "dq"
  )
  expect_error(
    cmBound(lambda.hat = 0.6385, p = 0.25, N = 310,
            dq = 0.073, N.dq = 1),
    "N.dq"
  )
})
