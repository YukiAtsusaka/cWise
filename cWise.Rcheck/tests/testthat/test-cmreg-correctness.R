test_that("outcome-model estimates recover known parameters", {
  set.seed(20260825)
  n <- 20000L
  x <- stats::rnorm(n)
  truth <- c(beta_intercept = -0.4, beta_x = 0.5,
             theta_intercept = 1.0, theta_x = -0.3)
  p <- 0.2
  p.prime <- 0.2
  design <- cbind(1, x)
  prevalence <- stats::plogis(drop(design %*% truth[1:2]))
  attention <- stats::plogis(drop(design %*% truth[3:4]))
  y_probability <- 0.5 + ((2 * p - 1) * prevalence + 0.5 - p) * attention
  a_probability <- 0.5 + (0.5 - p.prime) * attention
  simulated <- data.frame(
    Y = stats::rbinom(n, 1, y_probability),
    x = x,
    A = stats::rbinom(n, 1, a_probability)
  )

  fit <- cmreg(Y ~ x, anchor = A, p = p, p.prime = p.prime,
               data = simulated, n.start = 1L, control = list(maxit = 1200))
  expect_lt(max(abs((fit$estimates - truth) / fit$std.errors)), 3)
})

test_that("known latent-trait Gaussian block agrees with OLS", {
  set.seed(20260826)
  n <- 1000L
  x <- stats::rnorm(n)
  z <- stats::rbinom(n, 1, 0.35)
  outcome <- 0.5 + 0.8 * x + 1.2 * z + stats::rnorm(n, sd = 0.7)
  ols <- stats::lm(outcome ~ x + z)
  normal_loglik <- function(par) {
    mean <- par[1] + par[2] * x + par[3] * z
    sum(stats::dnorm(outcome, mean = mean, sd = exp(par[4]), log = TRUE))
  }
  reduced_fit <- stats::optim(c(0, 0, 0, 0), normal_loglik,
                               method = "BFGS", control = list(fnscale = -1))

  expect_equal(unname(reduced_fit$par[1:3]), unname(stats::coef(ols)), tolerance = 1e-6)
  expect_equal(exp(reduced_fit$par[4]), sqrt(mean(stats::residuals(ols)^2)), tolerance = 1e-6)
})

test_that("analytic outcome Hessian matches finite differences", {
  set.seed(20260827)
  n <- 200L
  x <- stats::rnorm(n)
  data <- data.frame(
    Y = stats::rbinom(n, 1, stats::plogis(-0.5 + 0.4 * x)),
    x = x,
    A = stats::rbinom(n, 1, stats::plogis(0.8 - 0.2 * x))
  )
  fit <- cmreg(Y ~ x, anchor = A, p = 0.2, p.prime = 0.2,
               data = data, n.start = 1L)
  model_frame <- stats::model.frame(Y ~ x, data)
  design <- stats::model.matrix(Y ~ x, model_frame)
  loglik <- function(par) {
    cWise:::cm_outcome_loglik(par, design, data$Y, data$A, 0.2, 0.2)
  }
  finite_difference_hessian <- function(par, step = 1e-4) {
    npar <- length(par)
    result <- matrix(0, npar, npar)
    for (i in seq_len(npar)) for (j in seq_len(npar)) {
      ei <- ej <- rep(0, npar)
      ei[i] <- step
      ej[j] <- step
      result[i, j] <- (loglik(par + ei + ej) - loglik(par + ei - ej) -
        loglik(par - ei + ej) + loglik(par - ei - ej)) / (4 * step^2)
    }
    result
  }
  numeric_hessian <- finite_difference_hessian(fit$estimates)
  relative_error <- max(abs(numeric_hessian - fit$hessian) /
    (1 + abs(fit$hessian)))
  expect_lt(relative_error, 1e-4)
})

test_that("shipped-data regression fixture remains stable", {
  outcome_fit <- cmreg(Y ~ female + age, anchor = A, p = 0.1, p.prime = 0.15,
                       data = cmdata2, n.start = 1L)
  predictor_fit <- cmreg_p(V ~ age + female, crosswise = Y, anchor = A,
                           p = 0.1, p.prime = 0.15, data = cmdata3,
                           n.start = 1L)
  expect_equal(
    unname(outcome_fit$estimates),
    c(-1.6501039, 0.2812883, 0.0326249, 0.1405443, -0.2058741, 0.0594322),
    tolerance = 1e-4
  )
  # The sigma element and predictor intervals intentionally changed during the
  # Phase 2 index-map and log-sigma corrections.
  expect_equal(
    unname(predictor_fit$estimates),
    c(-1.7342403, 0.0351695, 0.2879769, 0.2468424, 0.0548338,
      -0.1076469, 0.0235680, 0.0096394, 0.2472768, 0.9859290, 1.0438144),
    tolerance = 1e-4
  )
})

test_that("formula-general outcome and predictor fits work", {
  set.seed(20260828)
  outcome_data <- cmdata2
  outcome_data$group <- factor(sample(c("a", "b", "c"), nrow(outcome_data), TRUE))
  outcome_data$age[1] <- NA
  expect_s3_class(cmreg(Y ~ female, anchor = A, p = .1, p.prime = .15,
                         data = outcome_data, n.start = 1L), "cmreg")
  expect_s3_class(cmreg(Y ~ 1, anchor = A, p = .1, p.prime = .15,
                         data = outcome_data, n.start = 1L), "cmreg")
  expect_s3_class(cmreg(Y ~ group + age, anchor = A, p = .1, p.prime = .15,
                         data = outcome_data, n.start = 1L), "cmreg")
  expect_s3_class(cmreg(Y ~ female * age, anchor = A, p = .1, p.prime = .15,
                         data = outcome_data, n.start = 1L), "cmreg")
  expect_s3_class(cmreg(Y ~ poly(age, 2), anchor = A, p = .1, p.prime = .15,
                         data = outcome_data, n.start = 1L), "cmreg")

  predictor_data <- cmdata3
  predictor_data$group <- factor(sample(c("a", "b", "c"), nrow(predictor_data), TRUE))
  expect_s3_class(cmreg_p(V ~ group + age, crosswise = Y, anchor = A,
                           p = .1, p.prime = .15, data = predictor_data,
                           n.start = 1L), "cmreg_p")
})

test_that("degenerate inputs fail clearly", {
  expect_error(cmreg(Y ~ female, anchor = A, p = .5, p.prime = .15,
                     data = cmdata2), "must not equal 0.5")
  expect_error(cmreg(Y ~ female, anchor = A, p = 1.1, p.prime = .15,
                     data = cmdata2), "strictly between")
  constant <- cmdata2
  constant$Y <- 0
  expect_error(cmreg(Y ~ female, anchor = A, p = .1, p.prime = .15,
                     data = constant), "must not be constant")
  tiny <- data.frame(Y = c(0, 1, 0, 1), x = c(0, 1, 0, 1), A = c(0, 1, 1, 0))
  expect_error(cmreg(Y ~ x, anchor = A, p = .1, p.prime = .15, data = tiny),
               "must exceed")
  separated <- data.frame(Y = rep(c(0, 1), each = 10), x = rep(c(0, 1), each = 10),
                          A = rep(c(0, 1), 10))
  expect_error(cmreg(Y ~ x, anchor = A, p = .1, p.prime = .15, data = separated),
               "perfectly separated")
  constant_outcome <- cmdata3
  constant_outcome$V <- 1
  expect_error(cmreg_p(V ~ age + female, crosswise = Y, anchor = A,
                       p = .1, p.prime = .15, data = constant_outcome),
               "must not be constant")
})
