# Central parameter layout for the crosswise regression models.
#
# Keep all parameter slicing in one place. The predictor model estimates
# log_sigma internally; `sigma` identifies that same position on the reporting
# scale after the delta transformation.
cm_par_index <- function(k, model = c("outcome", "predictor")) {
  model <- match.arg(model)

  if (length(k) != 1L || is.na(k) || k < 1L || k != as.integer(k)) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }
  k <- as.integer(k)

  beta <- seq_len(k)
  theta <- k + seq_len(k)

  if (model == "outcome") {
    return(list(
      beta = beta,
      theta = theta,
      gamma = integer(),
      gamma_z = integer(),
      sigma = integer(),
      log_sigma = integer(),
      npar = 2L * k
    ))
  }

  gamma <- 2L * k + seq_len(k)
  gamma_z <- 3L * k + 1L
  sigma <- gamma_z + 1L

  list(
    beta = beta,
    theta = theta,
    gamma = gamma,
    gamma_z = gamma_z,
    sigma = sigma,
    # Reserved alias for the log-scale reparameterization in task 2.3.
    log_sigma = sigma,
    npar = sigma
  )
}

# Evaluate a column argument supplied either bare (for example, anchor = A)
# or as a single column name. Keeping this separate from the formula makes the
# anchor and crosswise responses explicit parts of the public API.
cm_data_column <- function(data, expression, argument, environment) {
  if (is.character(expression) && length(expression) == 1L) {
    if (!expression %in% names(data)) {
      stop(sprintf("`%s` must name a column in `data`.", argument), call. = FALSE)
    }
    value <- data[[expression]]
    name <- expression
  } else {
    value <- eval(expression, envir = data, enclos = environment)
    name <- paste(deparse(expression), collapse = "")
  }

  if (is.data.frame(value) || is.matrix(value)) {
    if (ncol(value) != 1L) {
      stop(sprintf("`%s` must evaluate to one column.", argument), call. = FALSE)
    }
    value <- value[, 1L]
  }
  if (length(value) != nrow(data)) {
    stop(sprintf("`%s` must have one value per row of `data`.", argument), call. = FALSE)
  }

  list(value = value, name = name)
}

cm_glm_start <- function(x, y) {
  fit <- tryCatch(
    stats::glm.fit(x = x, y = y, family = stats::binomial()),
    error = function(e) NULL
  )
  coefficients <- if (is.null(fit)) rep(0, ncol(x)) else fit$coefficients
  coefficients[!is.finite(coefficients)] <- 0
  pmax(pmin(coefficients, 5), -5)
}

cm_lm_start <- function(x, y) {
  fit <- tryCatch(stats::lm.fit(x = x, y = y), error = function(e) NULL)
  coefficients <- if (is.null(fit)) rep(0, ncol(x)) else fit$coefficients
  coefficients[!is.finite(coefficients)] <- 0
  coefficients
}

cm_hessian_vcov <- function(hessian, npar) {
  unavailable <- function(reason) {
    warning(
      sprintf("Standard errors are unavailable because %s.", reason),
      call. = FALSE
    )
    matrix(NA_real_, nrow = npar, ncol = npar)
  }

  if (!is.matrix(hessian) || !identical(dim(hessian), c(npar, npar)) ||
      any(!is.finite(hessian))) {
    return(unavailable("the Hessian is missing, malformed, or non-finite"))
  }

  hessian <- (hessian + t(hessian)) / 2
  eigenvalues <- tryCatch(
    eigen(hessian, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NULL
  )
  if (is.null(eigenvalues) || any(!is.finite(eigenvalues)) ||
      any(eigenvalues >= -sqrt(.Machine$double.eps))) {
    return(unavailable("the Hessian is not negative definite"))
  }

  inverse <- tryCatch(solve(hessian), error = function(e) NULL)
  if (is.null(inverse) || any(!is.finite(inverse))) {
    return(unavailable("the Hessian cannot be inverted"))
  }

  -inverse
}

cm_multistart_optim <- function(fn, start, n.start, control) {
  if (length(n.start) != 1L || is.na(n.start) || n.start < 1L ||
      n.start != as.integer(n.start)) {
    stop("`n.start` must be a single positive integer.", call. = FALSE)
  }
  if (!is.list(control)) {
    stop("`control` must be a list passed to `optim()`.", call. = FALSE)
  }

  n.start <- as.integer(n.start)
  optim_control <- utils::modifyList(list(maxit = 800, fnscale = -1), control)
  if (length(optim_control$fnscale) != 1L ||
      !isTRUE(all.equal(optim_control$fnscale, -1))) {
    warning("`control$fnscale` is ignored; cWise maximizes its log-likelihood.",
            call. = FALSE)
    optim_control$fnscale <- -1
  }

  starts <- vector("list", n.start)
  starts[[1L]] <- start
  if (n.start > 1L) {
    for (i in 2:n.start) {
      starts[[i]] <- start + stats::rnorm(length(start), sd = 0.1)
    }
  }

  fits <- lapply(starts, function(candidate) {
    tryCatch(
      stats::optim(
        par = candidate, fn = fn, method = "BFGS", control = optim_control,
        hessian = TRUE
      ),
      error = function(e) e
    )
  })
  valid <- vapply(fits, function(fit) {
    inherits(fit, "list") && is.finite(fit$value) && all(is.finite(fit$par))
  }, logical(1))
  converged <- valid & vapply(fits, function(fit) {
    inherits(fit, "list") && identical(fit$convergence, 0L)
  }, logical(1))

  if (!any(converged)) {
    stop(
      sprintf("Optimization failed to converge from all %d start(s).", n.start),
      call. = FALSE
    )
  }
  failed_starts <- n.start - sum(converged)
  if (failed_starts > 0L) {
    warning(
      sprintf("%d of %d optimization start(s) failed or did not converge; using the best converged fit.",
              failed_starts, n.start),
      call. = FALSE
    )
  }

  candidate_indices <- which(converged)
  fits[[candidate_indices[which.max(vapply(fits[candidate_indices], `[[`, numeric(1), "value"))]]]
}

cm_prediction_data <- function(out, newdata, zval, typical) {
  if (!is.null(newdata) && !is.null(typical)) {
    stop("Supply either `newdata` or `typical`, not both.", call. = FALSE)
  }
  supplied <- if (is.null(newdata)) typical else newdata
  argument <- if (is.null(newdata)) "typical" else "newdata"

  if (is.null(supplied)) {
    values <- list()
  } else if (is.data.frame(supplied)) {
    values <- as.list(supplied)
  } else if (is.list(supplied) && !is.null(names(supplied))) {
    values <- supplied
  } else if (is.atomic(supplied) && !is.null(names(supplied))) {
    values <- as.list(supplied)
  } else {
    stop(sprintf("`%s` must be a named data frame, list, or vector.", argument),
         call. = FALSE)
  }

  rhs_terms <- stats::delete.response(out$terms)
  required <- all.vars(rhs_terms)
  if (!is.null(zval)) {
    z_names <- unique(names(zval))
    if (is.null(names(zval)) || length(z_names) != 1L || !nzchar(z_names)) {
      stop("`zval` must be a named vector for exactly one covariate.", call. = FALSE)
    }
    values[[z_names]] <- unname(zval)
  }

  unknown <- setdiff(names(values), required)
  if (length(unknown)) {
    stop(sprintf("Unknown prediction variable(s): %s.",
                 paste(unknown, collapse = ", ")), call. = FALSE)
  }
  missing <- setdiff(required, names(values))
  if (length(missing)) {
    stop(sprintf("Missing prediction variable(s): %s.",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }

  if (!length(required)) {
    prediction_data <- data.frame(.cm_row = 1L)
  } else {
    lengths <- vapply(values[required], length, integer(1))
    if (any(lengths == 0L)) {
      stop("Prediction variables cannot be empty.", call. = FALSE)
    }
    n <- max(lengths)
    if (any(!lengths %in% c(1L, n))) {
      stop("Prediction variables must have length one or a common length.",
           call. = FALSE)
    }
    values <- lapply(values[required], rep, length.out = n)
    prediction_data <- as.data.frame(values, optional = TRUE)
  }

  for (variable in intersect(names(out$xlevels), names(prediction_data))) {
    if (!is.factor(prediction_data[[variable]])) {
      prediction_data[[variable]] <- factor(
        prediction_data[[variable]], levels = out$xlevels[[variable]]
      )
    }
  }

  design <- tryCatch(
    stats::model.matrix(
      rhs_terms, data = prediction_data, xlev = out$xlevels,
      contrasts.arg = out$contrasts
    ),
    error = function(e) {
      stop("Could not construct the prediction design matrix: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  missing_columns <- setdiff(out$design_columns, colnames(design))
  if (length(missing_columns)) {
    stop(sprintf("Prediction design is missing column(s): %s.",
                 paste(missing_columns, collapse = ", ")), call. = FALSE)
  }
  design <- design[, out$design_columns, drop = FALSE]
  if (any(!is.finite(design))) {
    stop("`newdata` contains missing or non-finite values.", call. = FALSE)
  }
  design
}

cm_with_seed <- function(seed, generator) {
  if (is.null(seed)) {
    return(generator())
  }
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed != as.integer(seed)) {
    stop("`seed` must be a single integer or `NULL`.", call. = FALSE)
  }

  seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (seed_exists) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (seed_exists) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  generator()
}

cm_prediction_summary <- function(draws, estimate, include_draws) {
  result <- data.frame(
    estimate = as.numeric(estimate),
    conf.low = apply(draws, 1L, stats::quantile, probs = 0.025, names = FALSE),
    conf.high = apply(draws, 1L, stats::quantile, probs = 0.975, names = FALSE)
  )
  rownames(result) <- rownames(draws)
  if (isTRUE(include_draws)) {
    attr(result, "draws") <- draws
  }
  result
}
