## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----setup--------------------------------------------------------------------
library(cWise)

## ----trait-outcome------------------------------------------------------------
outcome_fit <- cmreg(
  Y ~ female + age,
  anchor = A, p = 0.10, p.prime = 0.15,
  data = cmdata2, n.start = 1L
)
summary(outcome_fit)

## ----trait-prediction---------------------------------------------------------
trait_predictions <- cmpredict(
  outcome_fit,
  newdata = data.frame(female = c(0, 1), age = c(30, 30)),
  nsim = 300, seed = 20260825
)
trait_predictions

## ----trait-prediction-plot, fig.width=6, fig.height=4-------------------------
trait_labels <- c("Female = 0", "Female = 1")
plot(
  seq_len(nrow(trait_predictions)), trait_predictions$estimate,
  ylim = c(0, 1), xlim = c(0.75, 2.25), xaxt = "n",
  xlab = "Covariate scenario", ylab = "Predicted trait prevalence",
  pch = 19, col = "#1B6CA8"
)
axis(1, at = seq_along(trait_labels), labels = trait_labels)
arrows(
  x0 = seq_len(nrow(trait_predictions)), y0 = trait_predictions$conf.low,
  x1 = seq_len(nrow(trait_predictions)), y1 = trait_predictions$conf.high,
  angle = 90, code = 3, length = 0.05, col = "#1B6CA8"
)

## ----trait-predictor----------------------------------------------------------
predictor_fit <- cmreg_p(
  V ~ age + female,
  crosswise = Y, anchor = A, p = 0.10, p.prime = 0.15,
  data = cmdata3, n.start = 1L
)
summary(predictor_fit)

## ----outcome-prediction-------------------------------------------------------
outcome_predictions <- cmpredict_p(
  predictor_fit,
  newdata = data.frame(age = 30, female = 1),
  nsim = 300, seed = 20260825
)
outcome_predictions

## ----outcome-prediction-plot, fig.width=6, fig.height=4-----------------------
state_labels <- sub("^[0-9]+: ", "", rownames(outcome_predictions))
ylim <- range(outcome_predictions$conf.low, outcome_predictions$conf.high)
plot(
  seq_len(nrow(outcome_predictions)), outcome_predictions$estimate,
  ylim = ylim + c(-0.05, 0.05), xlim = c(0.75, 2.25), xaxt = "n",
  xlab = "Latent-trait scenario", ylab = "Predicted observed outcome",
  pch = 19, col = "#A23A2E"
)
axis(1, at = seq_along(state_labels), labels = state_labels)
arrows(
  x0 = seq_len(nrow(outcome_predictions)), y0 = outcome_predictions$conf.low,
  x1 = seq_len(nrow(outcome_predictions)), y1 = outcome_predictions$conf.high,
  angle = 90, code = 3, length = 0.05, col = "#A23A2E"
)

