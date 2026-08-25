# Central parameter layout for the crosswise regression models.
#
# Keep all parameter slicing in one place. The predictor model still estimates
# sigma on its natural scale; task 2.3 will change that parameter to log_sigma
# without changing its position in the vector.
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
