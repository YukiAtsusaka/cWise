# cWise — CRAN Submission Plan

Working document for **Yuki Atsusaka** and **Kolbe Dumas**. Update the status boxes as
you go. The Phase 0 check output was verified against the package as of commit `a27f64f`
(R 4.4.1, Windows) — the check output in Phase 0 is real, not hypothetical.

**Target:** `cWise 0.1.0` on CRAN, with two vignettes.

**Progress through 2026-08-25:**

- Completed 1.1: migrated `DESCRIPTION` to `Authors@R`; Yuki Atsusaka is `aut`/`cre`,
  and Kolbe Dumas and Randy T. Stevenson are `aut`. Kolbe's package email is
  `kdumas@uh.edu`.
- Maintainer decision: YA is the package maintainer, represented by the `cre` role in
  `Authors@R`; R derives the `Maintainer` field from Yuki's name and email.
- Completed 1.15: replaced the dead working-paper citation with the published 2023
  *Political Analysis* article and DOI `10.1017/pan.2021.43`.
- Completed 1.3: the title is now "Crosswise Models for Sensitive Survey Questions";
  the description is two sentences and cites the DOI in CRAN format.
- Completed 1.4: bumped the development version to `0.1.0` and rewrote `NEWS.md`
  to record the completed CRAN-submission preparation work.
- Completed 1.5: verified that the standard `License: GPL-3` declaration is valid;
  a separate `LICENSE` file is neither needed nor appropriate unless the declaration
  explicitly refers to one.
- Completed 1.6: declared the missing `stats`, `graphics`, and `utils` functions in
  `R/cWise-package.R` and regenerated `NAMESPACE` and internal package help with
  roxygen 7.3.3. A fresh RStudio tarball check reported zero missing declarations.
- Completed 1.8: `bc.est()` now uses explicit `digits =` arguments, and no
  abbreviated `round(..., d = ...)` calls remain after the regression rewrite.
- Completed 1.9: removed the legacy `R/power_analysis.R` file, whose duplicate
  power functions were silently superseded and whose undocumented `sim.curve()`
  helper was broken. A fresh RStudio build and install retained `sim.power()` and
  `sim.power.N()` and no longer contained `sim.curve()`.
- Updated the README and simulation references from the pre-publication year to 2023.
- Added the two-person ownership map under "Suggested sequencing for two people."

**Owner column:** put initials in the `Who` cell when you pick a task up. `YA` = Yuki,
`KD` = Kolbe, `?` = unassigned.

---

## How to run the checks (do this before touching anything)

```r
# from the package root
devtools::document()                      # regenerate NAMESPACE + man/ from roxygen
devtools::load_all()                      # quick interactive reload
devtools::test()                          # testthat
devtools::check(args = "--as-cran")       # the gate that matters
```

Two extra gates required before submission:

```r
devtools::check_win_devel()   # emails you the Windows r-devel result
rhub::rhub_check()            # Linux + macOS; needs a one-time rhub::rhub_setup()
```

Rule of thumb: **CRAN accepts 0 ERRORs, 0 WARNINGs, and only NOTEs you can justify
in the submission comments.** The only NOTE we expect to keep is `New submission`.

---

## Phase 0 — Baseline: what `R CMD check --as-cran` says today

This is the original baseline snapshot, not a live check report. The check aborted
grows as we fix things — expect new findings after each phase.

**ERROR — `DESCRIPTION`**
```
Required field missing or empty: 'Maintainer'
```
**Resolved 2026-08-25:** the invalid `Author: person(...)` field was replaced with a
valid `Authors@R` list. R now derives Yuki Atsusaka as maintainer from the `cre` role.
Re-run the full check after the remaining Phase 1 metadata changes to confirm this
gate stays resolved.

**ERROR — dependencies** (appears once `DESCRIPTION` parses)
```
VignetteBuilder package not declared: 'knitr'
Namespace dependencies missing from DESCRIPTION Imports/Depends entries: 'dplyr', 'ggplot2'
```
`NAMESPACE` imports from `dplyr` and `ggplot2`, but `DESCRIPTION` lists `tidyverse`.
**CRAN will not accept `tidyverse` as an `Imports` entry** — it is a meta-package.

**Resolved 2026-08-25:** `Imports` now declares `dplyr`, `ggplot2`, `mvtnorm`, and
`scales` individually; `knitr` and `rmarkdown` are in `Suggests`, and unused `rlang`
was removed. A source-package check reached `checking package dependencies ... OK`.

**ERROR — examples**
```
Error in ggtitle("Sensitivity Analysis") : could not find function "ggtitle"
```
`man/cmBound.Rd` calls `ggtitle()`/`theme()`/`element_text()` in the example, but
ggplot2 is imported, not attached, so those names are not visible to the user's session.
The check dies here, which means **the `cmpredict`/`cmreg` examples have not run yet**
and their failures are still hidden. Known bad ones:
- `man/cmpredict.Rd`: `m <- cmpredict(m, zval=0, typical=30)` — `m` is undefined.
- `man/cmpredict.p.Rd`: `pr2 <- cmpredict2(out=m2, typical=c(1,30))` — `cmpredict2()`
  does not exist, and `m2` is undefined.

**ERROR — tests**
```
Error in test_dir(): No test files found
```
`tests/testthat.R` exists but `tests/testthat/` does not.

**WARNING — install-time**
```
Note: possible error in 'sim.power(N.sim = N.sim, ...)': argument 3 matches multiple formal arguments
```
Root cause: **`sim.power()` and `sim.power.N()` are each defined twice**, with
conflicting signatures.

| Function | `R/power_analysis.R` | `R/sim_power.R` / `R/sim_power_N.R` |
|---|---|---|
| `sim.power` | `(N.sim, sample, pi, p, p.prime, gamma, direct)` | `(N.sim, sample, pi.null, pi.alt, p, p.prime, gamma, direct)` |
| `sim.power.N` | `(N.sim, pi, p, p.prime, gamma, direct)` | `(N.sim, pi, p, p.prime, gamma, direct)` |

Files load alphabetically, so `sim_power.R` wins and `power_analysis.R`'s versions are
silently discarded. `sim.curve()` (only defined in `power_analysis.R`, exported nowhere,
documented nowhere) then calls `sim.power(pi = ...)`, which partially matches *both*
`pi.null` and `pi.alt` → runtime error. **`sim.curve()` is currently dead and broken.**

**WARNING — undeclared / unused imports**
```
'::' or ':::' import not declared from: 'scales'
Namespace in Imports field not imported from: 'rlang'
```

**WARNING — `\usage` mismatch**
```
Documented arguments not in \usage in Rd file 'cmpredict.p.Rd': 'zval'
```
The roxygen block documents `zval`, but `cmpredict.p()` takes `(out, typical)`.

**NOTEs**
- `Version contains large components (0.0.0.9000)` — dev version; ship `0.1.0`.
- ~~`inst/CITATION` URL returns **404**.~~ Resolved 2026-08-25 by replacing the
  working-paper entry with the published article and DOI.
- ~45 × `no visible global function definition` — every base-package function
  (`optim`, `pnorm`, `quantile`, `model.frame`, `rbinom`, `abline`, `txtProgressBar`, …)
  needs an `importFrom`. The check prints the exact block to paste; see Phase 1.
- `partial argument match of 'd' to 'digits'` in `round(x, d = 4)` throughout.

---

## Phase 1 — Mechanical CRAN compliance

Cheap, independent, no statistics involved. Good place for Kolbe to start
while the estimation work (Phase 2) is in flight. Nothing here changes numerical output.

| # | Task | File(s) | Who | Done |
|---|---|---|---|---|
| 1.1 | Replace `Author:` with `Authors@R = c(person(..., role=c("aut","cre")), person("Randy","Stevenson", role="aut"))`. Drop the `Author:`/`Maintainer:` fields entirely. Decide whether Kolbe is listed — `role="aut"` if he contributes package code, `role="ctb"` for smaller contributions. | `DESCRIPTION` | KD | ☒ |
| 1.2 | Replace `tidyverse` in `Imports` with the packages actually used: `dplyr`, `ggplot2`, `scales`, `mvtnorm`. Drop `rlang` unless 1.9 keeps it. Add `knitr`, `rmarkdown` to `Suggests`. | `DESCRIPTION` | KD | ☒ |
| 1.3 | Rewrite `Title` in title case, no redundant article, ≤65 chars. `Description` must be ≥2 sentences, must not start with "This package", and must cite the paper as `Atsusaka and Stevenson (2023) <doi:...>`. | `DESCRIPTION` | KD | ☒ |
| 1.4 | Bump `Version: 0.1.0`. Rewrite `NEWS.md` under a `# cWise 0.1.0` heading (it currently says `0.0.1`, which matches nothing). | `DESCRIPTION`, `NEWS.md` | KD | ☒ |
| 1.5 | Verified standard `License: GPL-3`; no separate `LICENSE` file should be added unless the declaration refers to it. | `DESCRIPTION` | KD | ☒ |
| 1.6 | Declared the full `stats`, `graphics`, and `utils` import set in `R/cWise-package.R`; regenerated `NAMESPACE` and package help with `devtools::document()`. A fresh tarball check found zero missing declarations. | `R/cWise-package.R`, `NAMESPACE`, `man/cWise-package.Rd` | KD | ☒ |
| 1.7 | Add `@importFrom scales ...` or fully qualify the `scales::` calls, and declare `scales` in `Imports`. Calls in `R/sim_estimates.R` are fully qualified. | `R/sim_estimates.R`, `DESCRIPTION` | KD | ☒ |
| 1.8 | Confirmed that no abbreviated `round(..., d = ...)` calls remain after the regression rewrite; `bc.est()` now uses explicit `digits =` arguments. | `R/bc.est.R`, `R/cmreg.R`, `R/cmreg.p.R` | KD | ☒ |
| 1.9 | Removed legacy `R/power_analysis.R` after confirming that `sim_power.R` and `sim_power_N.R` provide the active supported functions. The undocumented `sim.curve()` helper was broken and removed; `NEWS.md` records the change. A fresh build and install passed. | `R/power_analysis.R`, `NEWS.md` | KD | ☒ |
| 1.10 | Completed as part of Phase 2.3; no global warning option is changed. | `R/cmreg.p.R` | YA | ☒ |
| 1.11 | Remove the hard-coded `set.seed()` calls from inside exported functions. They silently overwrite the user's RNG stream, which is a CRAN policy problem and makes results non-reproducible in the way users expect. Add a `seed = NULL` argument and, when set, restore state on exit via `on.exit()`. **Progress:** completed for `cmpredict()` and `cmpredict.p()` in 2.8; `bc.est()` remains. | `R/bc.est.R:63`, `R/cmpredict.R`, `R/cmpredict.p.R` | KD | ☐ |
| 1.12 | Route all `cat()`/`print()` progress output through `message()` and gate it behind a `verbose = TRUE` argument, so it can be silenced and goes to stderr. | `R/sim_power.R`, `R/sim_cwdata.R` | ? | ☐ |
| 1.13 | Fix the `ggtitle` example: either `\dontrun{}` the ggplot2-syntax lines or prefix `library(ggplot2)` inside the example. Prefer the latter — a runnable example is better. | `R/cmBound.R` roxygen | ? | ☐ |
| 1.14 | Replace `size` with `linewidth` in `geom_line`/`geom_segment` calls (deprecated since ggplot2 3.4.0; currently emits a deprecation warning during checks). | `R/cmBound.R` | ? | ☐ |
| 1.15 | Fix the `inst/CITATION` 404 URL and update `year`/`journal` to the published reference. Use `bibentry(bibtype = "Article", ...)` — capital A. | `inst/CITATION` | KD | ☒ |
| 1.16 | Add `^doc$`, `^\.github$`, `^README\.Rmd$`, `^cran-comments\.md$`, `^CLAUDE\.md$` to `.Rbuildignore`. Decide whether `doc/Appendix_C5.pdf` ships in `inst/` (it is a real reference, so probably yes → move it to `inst/doc_ref/`, **not** `inst/doc/`, which knitr owns). | `.Rbuildignore` | ? | ☐ |
| 1.17 | Consider renaming the `pi` argument in the `sim.*` functions to `prev` or `pi.true`. `pi` masks the base constant inside those function bodies and is what caused the 1.9 partial-matching bug. Not a check failure today — a trap for later. | `R/sim_*.R` | ? | ☐ |

**Exit criterion:** `devtools::check(args = "--as-cran")` reaches the examples/tests
stage and reports only NOTEs about `New submission` plus whatever Phase 2 breaks.

---

## Phase 2 — The estimation and prediction rewrite (the real blocker)

This is the substantive work, and it should not be rushed to hit a submission date.
The functions currently work *only* for the shapes used in the shipped examples.

### 2a. What is actually wrong (all reproduced, not speculative)

**(1) `cmreg()` breaks on any covariate that is not a single numeric column.**
`cmreg()` mixes two different widths: `X1` comes from `model.matrix()` (which expands
factors into dummies) while `n.var` comes from `dim(df)` on the `model.frame` (which does
not). They agree only when every predictor is one column.

```r
d <- cmdata2; d$grp <- factor(sample(c("a","b","c"), nrow(d), TRUE))
cmreg(Y ~ female + grp + A, p = 0.1, p.prime = 0.15, data = d)
#> Error: undefined columns selected
```

So factors, interactions, `poly()`, and `I()` all fail. Same defect in `cmreg.p()`.

**(2) `cmreg.p()` mislabels its own output: the row printed as `sigma` is the gamma intercept.**
The residual SD is `par[k.g.e]` (index 11 in the shipped example), but `AuxiliaryCoef2`
is built from `par[(n.var+1):(2*n.var+1)]` = `par[4:7]`, and `par[7]` is `gamma0` — the
first element of the gamma block, which `Coefficients` also reports. Confirmed:

```r
m2 <- cmreg.p(V ~ age + female + Y + A, p = 0.1, p.prime = 0.15, data = cmdata3)
m2$AuxiliaryCoef2["sigma", 1]      #>  0.0235
m2$Coefficients["(intercept)", 1]  #>  0.0235   <- identical, same SE
```
**The true residual standard deviation is never reported anywhere.**

**(3) `cmpredict.p()` pairs a reordered mean vector with an unreordered covariance matrix.**
`optim()` estimates parameters in the order `(beta, theta, gamma, sigma)`, and `out$VCV`
preserves that order. But `cmpredict.p()` builds
`coefs <- c(coef.gamma, coef.beta, coef.theta)` and then calls
`rmvnorm(mean = coefs, sigma = vcovs)`. Element *i* of the mean no longer corresponds to
row/column *i* of the covariance matrix, so **the parametric-bootstrap draws are drawing
the wrong parameters against the wrong covariances.** This silently produces plausible-
looking but incorrect uncertainty intervals. Treat as the highest-severity item here.

**(4) Both prediction functions hard-code the design vector layout.**
- `cmpredict()`: `typ.vec <- cbind(1, zval, typical)` assumes the variable of interest is
  the *first* term in the formula and that `typical` is supplied in exactly the remaining
  order. Nothing validates the length; a wrong length either errors with a bare
  `non-conformable arguments` or, in the other direction, quietly returns a number.
- `cmpredict.p()`: `typ.vec <- rbind(c(1, typical, 0), c(1, typical, 1))` hard-codes
  exactly two scenarios and assumes the latent-trait coefficient is the last gamma element.
  It also returns the *linear predictor*, while `cmpredict()` returns a *probability* —
  inconsistent contracts for two functions documented identically.

**(5) The optimizer has no safety net.**
- `init <- rep(0.01, ...)` — one fixed start, not user-settable, no multi-start.
- `MLE$convergence` is never inspected. A non-converged fit is reported as if converged.
- `sigma` is estimated unconstrained by BFGS, so it can wander to ≤ 0 and put a negative
  value inside `log()`. This is almost certainly what `options(warn=-1)` was hiding.
- `SE <- sqrt(diag(-solve(H)))` assumes `H` is invertible and negative definite. If it is
  not, this yields `NaN` standard errors with no message, or `solve()` errors out.

**(6) Predictions are computed from *rounded* coefficients.**
`out$Coefficients[,1]` has already been through `round(digits = 4)`; `out$VCV` has not.
The fit object should carry full-precision estimates and round only at print time.

### 2b. Plan of attack

Do these in order. 2.1 is the foundation — it makes 2.2–2.4 straightforward instead of
another round of index arithmetic.

| # | Task | Who | Done |
|---|---|---|---|
| 2.1 | **Introduce an explicit parameter index map.** Added internal `cm_par_index()` and routed likelihoods, output assembly, VCV subsetting, and predictions through it. Removed literal parameter-range arithmetic. Verified unchanged `cmreg()` fixture estimates, corrected `cmreg.p()` sigma reporting (`1.0438` rather than gamma intercept `0.0235`), and aligned prediction means with matching VCV subsets. | YA | ☒ |
| 2.2 | **Completed 2026-08-25.** Widths now come exclusively from `model.matrix()`. `cmreg()` takes `anchor = A`; `cmreg.p()` takes `crosswise = Y, anchor = A`; these responses are no longer formula terms. Updated Rd and README examples and recorded the breaking change in `NEWS.md`. Verified factor, interaction, polynomial, intercept-only, and missing-value inputs for `cmreg()`, plus factor and missing-value inputs for `cmreg.p()`. | YA | ☒ |
| 2.3 | Reparameterized the predictor model as `log_sigma`: likelihood evaluation uses `sigma = exp(log_sigma)`, so the scale is always positive. Reported sigma, SE, and VCV are delta-transformed back to the sigma scale. Verified positive sigma (`1.0438`) and stable fixture estimates. | YA | ☒ |
| 2.4 | Added data-informed `start = NULL` defaults (binomial GLMs and least squares), reporting-scale user starts, `n.start = 3` multistart search, and `control` passthrough. The best converged likelihood is retained; partial failures warn and all-failed optimization stops. Hessians are checked for invertibility and negative definiteness; unavailable SEs/VCVs are reported as `NA` with a warning, never `NaN`. | YA | ☒ |
| 2.5 | Fit objects now retain full-precision `estimates`, `std.errors`, coefficient tables, named `VCV`, log-likelihood, `n`, convergence, `p`, `p.prime`, and the matched call. Added classes `"cmreg"` and `c("cmreg.p", "cmreg")`, plus registered S3 `print()`/`summary()` methods. Rounding occurs only during display. | YA | ☒ |
| 2.6 | `cmpredict()` now accepts named `newdata` or named `typical` inputs and optional named, vectorized `zval`. It rebuilds and aligns the design matrix from stored model-frame terms, factor levels, contrasts, and column names. Inputs are validated by name and length. Verified direct multi-row data, curves, factors, interactions, `poly()`, and intercept-only models. | YA | ☒ |
| 2.7 | `cmpredict.p()` now uses named `newdata`/`typical` inputs and vectorized named `zval`, returning absent/present latent-trait rows for each scenario. Bootstrap means and VCVs are subset through the same gamma-plus-gamma-z index set. Added `type = c("response", "link")`; they are equal under the current Gaussian identity-link model. Verified exact agreement with a manually indexed bootstrap draw, plus factor and invalid-input cases. | YA | ☒ |
| 2.8 | Both prediction functions now use `nsim = 1000`, `seed = NULL`, and `draws = FALSE` by default. Seeded calls restore the exact prior RNG state (including no prior seed). They return scenario-level `estimate`, `conf.low`, and `conf.high` data frames; `draws = TRUE` attaches the raw matrix. Verified reproducibility, interval calculations, and RNG preservation. | YA | ☒ |
| 2.9 | Renamed the canonical dotted API to `bc_est`, `cmreg_p`, `cmpredict_p`, `sim_cwdata`, `sim_estimates`, `sim_power`, and `sim_power_N`. Dotted names remain exported `.Deprecated()` wrappers for one release. The predictor fit class is now `c("cmreg_p", "cmreg")`; README, NEWS, examples, internal calls, docs, and NAMESPACE use the canonical names. | YA | ☒ |

### 2c. Correctness gates for Phase 2 — do not skip

The whole point of Phase 2 is that we cannot currently tell a correct fit from an
incorrect one. Add the evidence before declaring it fixed.

| # | Task | Who | Done |
|---|---|---|---|
| 2.10 | **Recovery test.** Simulate from the model with known parameters at large `n` (e.g. 20,000) and assert every estimate is within ~3 SEs of truth. This is the test that would have caught bugs (2) and (3). | YA | ☑ — `test-cmreg-correctness.R` |
| 2.11 | **Analytic cross-check on `sigma`.** With the latent trait held fixed, `cmreg.p`'s gamma/sigma block reduces to OLS; assert agreement with `lm()` to a tolerance. Pins down the residual SD that is currently unreported. | YA | ☑ — fixed-latent OLS comparison |
| 2.12 | **Numeric-vs-analytic Hessian check** on a small dataset, to confirm the SEs come from the likelihood we think we wrote. | YA | ☑ — central-difference comparison |
| 2.13 | **Regression fixture.** Freeze the current `cmdata2`/`cmdata3` coefficient estimates in a test file *before* refactoring, so 2.1–2.9 can be verified as behaviour-preserving except where we intend a change. Where output legitimately changes (the `sigma` row, `cmpredict.p` intervals), record the old and new values side by side in the commit message — those two are *supposed* to move. | YA | ☑ — canonical-fit coefficients frozen |
| 2.14 | **Formula-generality suite.** Fit with: one covariate; no covariates (intercept only — check `X0[, -ncol(X0)]` does not drop to a vector); a 3-level factor; an interaction; a `poly()` term; a covariate with `NA`s. All must fit or fail with a clear message. | YA | ☑ — includes transformed-formula NA handling |
| 2.15 | **Degenerate-input suite.** `p = 0.5` (division by `2*p-1` → `Inf`), `p` outside `(0,1)`, perfectly separated data, `n` smaller than the parameter count, constant `Y`. Each should `stop()` with an informative message rather than return nonsense. | YA | ☑ — validation/error tests |

---

## Phase 3 — Tests

`tests/testthat/` does not exist yet, so the suite currently ERRORs. Build it out
alongside Phase 2 rather than after.

| # | Task | Who | Done |
|---|---|---|---|
| 3.1 | `usethis::use_testthat(3)`; confirm `DESCRIPTION` gains `Config/testthat/edition: 3`. | ? | ☐ |
| 3.2 | `test-bc-est.R`: known-input snapshot; weighted vs unweighted; `NA` handling; that `weight` missing and `weight = rep(1,n)` agree. | ? | ☐ |
| 3.3 | `test-cmreg.R`, `test-cmreg-p.R`: the Phase 2c gates. | ? | ☐ |
| 3.4 | `test-cmpredict.R`: named-argument matching; informative error on wrong `typical` length; that predictions lie in `[0,1]`; that `seed` makes results reproducible and *does not* perturb the caller's RNG stream. | ? | ☐ |
| 3.5 | `test-cmbound.R`: returns a `ggplot` object; `kappa`/`dq`/`N.dq` variants; `p = 0.5` guarded. | ? | ☐ |
| 3.6 | `test-sim.R`: smoke tests at tiny `N.sim` so they run in seconds. Wrap anything slow in `skip_on_cran()`. | ? | ☐ |
| 3.7 | Keep total check time well under CRAN's 10-minute-per-platform budget. The `sim.*` functions are the risk — `skip_on_cran()` liberally and use `\donttest{}` in their examples. | ? | ☐ |

---

## Phase 4 — Vignettes

`vignettes/cWise_UserGuide.Rmd` is a stub: the only code chunk is `#library(cWise)`,
commented out. It builds, but it documents nothing.

Recommended split into two vignettes, because "estimate a prevalence" and "design a
study" are different audiences arriving with different questions:

| # | Task | Who | Done |
|---|---|---|---|
| 4.1 | `usethis::use_vignette("cWise")` — **"Getting started with cWise"**. Delete the stub. Narrative: what the crosswise model is; why inattentive respondents bias the naive estimator; `bc.est()` on `cmdata`; reading the output; `cmBound()` for sensitivity analysis. Should stand alone for a reader who has not read the paper. | ? | ☐ |
| 4.2 | `usethis::use_vignette("regression")` — **"Regression with a latent sensitive trait"**. `cmreg()` with the trait as outcome, `cmreg.p()` with it as predictor, then `cmpredict()`/`cmpredict.p()` for predicted probabilities with uncertainty. **Write this last** — it is the surface Phase 2 changes most, and rewriting it twice is wasted effort. | ? | ☐ |
| 4.3 | Decide whether power analysis gets a third vignette or a section in 4.1. If a vignette: pre-compute the simulations, save results to `inst/extdata/`, and load them in the chunk. A vignette that runs `sim.power(N.sim = 500)` at build time will blow the CRAN time budget. | ? | ☐ |
| 4.4 | Every vignette: `eval` must succeed from a clean session with only declared deps attached. No `tidyverse`, no internet access, no `setwd()`, no writing outside `tempdir()`. | ? | ☐ |
| 4.5 | Add `%\VignetteDepends{}` where needed; confirm `VignetteBuilder: knitr` and that `knitr`/`rmarkdown` are in `Suggests`. | ? | ☐ |
| 4.6 | `devtools::build_vignettes()`, then confirm `inst/doc` is gitignored (it already is) but *is* included in the built tarball. | ? | ☐ |
| 4.7 | Optional: convert `README.md` → `README.Rmd` so its output stays in sync with the code, and add `^README\.Rmd$` to `.Rbuildignore`. The current README is 12 KB of hand-maintained output. | ? | ☐ |

---

## Phase 5 — Docs polish

| # | Task | Who | Done |
|---|---|---|---|
| 5.1 | Every exported function needs `@return` describing the object concretely (`\describe{}` over the list elements), plus `@examples` that **actually run**. Fix `cmpredict`/`cmpredict.p` examples so they build `m`/`m2` first. | ? | ☐ |
| 5.2 | Reconcile the `cmpredict.p` roxygen block with its real signature (drop `zval` or add the argument). | ? | ☐ |
| 5.3 | Add `R/cWise-package.R` with `@keywords internal` + `"_PACKAGE"` for the package-level `?cWise` page, and hang the base `@importFrom` tags there. | ? | ☐ |
| 5.4 | Cross-link with `@seealso`; add `@references` with the paper DOI once available. | ? | ☐ |
| 5.5 | Spell check: `devtools::spell_check()`. | ? | ☐ |
| 5.6 | Optional but useful: `usethis::use_github_action("check-standard")` so both of us see check results on every push instead of discovering breakage at submission time. | ? | ☐ |

---

## Phase 6 — Submission

| # | Task | Who | Done |
|---|---|---|---|
| 6.1 | `devtools::check(args = "--as-cran")` clean locally. | ? | ☐ |
| 6.2 | `devtools::check_win_devel()` — clean on Windows r-devel. | ? | ☐ |
| 6.3 | `rhub::rhub_check()` — clean on Linux and macOS. | ? | ☐ |
| 6.4 | Write `cran-comments.md`: platforms tested, R versions, and a one-line justification for each remaining NOTE (expect only `New submission`). Add it to `.Rbuildignore`. | ? | ☐ |
| 6.5 | Confirm the maintainer email is one Yuki will hold long-term and can respond from within a few days — CRAN archives packages whose maintainer goes unreachable. YA is assigned as maintainer; confirm the current `atsusaka@uh.edu` address is durable before submission. | YA | ☐ |
| 6.6 | Tag the release in git (`v0.1.0`) and confirm the working tree is clean. | ? | ☐ |
| 6.7 | `devtools::release()` — walks the final checklist and submits. Yuki must confirm from the maintainer address. | ? | ☐ |
| 6.8 | Expect a human reviewer round. Reply on the same email thread, bump to `0.1.1`, resubmit. Budget 1–3 weeks. | ? | ☐ |

---

## Suggested sequencing for two people

Phase 1 and Phase 2 are almost fully independent — one is metadata and namespace
hygiene, the other is the likelihood code. Splitting them is the main parallelism win.

1. **Both:** agree on the Phase 2.2 and 2.9 API decisions first. They are breaking
   changes, and everything else (docs, tests, vignettes) is written against them.
   Doing this out of order is the most likely source of duplicated work.
2. **Kolbe:** Phase 1 (all mechanical), then Phase 4.1 (getting-started vignette —
   it only touches `bc.est`/`cmBound`, which Phase 2 does not change), then Phase 5.
3. **Yuki:** Phase 2 (estimation rewrite) with Phase 2c gates and Phase 3.3–3.4 tests
   written as you go.
4. **Rejoin:** Phase 4.2 (regression vignette) once the Phase 2 API is frozen.
5. **Both:** Phase 6.

Do not start Phase 6 until Phase 2c is signed off. Submitting with bug (3) — the
mean/covariance mismatch in `cmpredict.p()` — would publish incorrect confidence
intervals under our names.

---

### Assigned work

Use this as the ownership map for the `Who` column above. Shared items require both
people to review and sign off; the first person listed is the implementation lead.

| Owner | Assigned tasks | Dependency / handoff |
|---|---|---|
| YA + KD | Settle open question 3 | Decisions 2.2 and 2.9 and the paper citation were approved/resolved on 2026-08-25. Release scope remains open. |
| KD | 1.1-1.17 | Start immediately. Ask YA only for author/maintainer, simulation-scope, and citation decisions. |
| YA | 2.1-2.9 | Start after the joint API decisions. Keep the index-map and full-precision fit-object work ahead of prediction rewrites. |
| YA | 2.10-2.15, 3.3-3.4 | Write these correctness and regression tests alongside the corresponding Phase 2 changes, not afterward. |
| KD | 3.1-3.2, 3.5-3.7 | Establish testthat first, then cover the mechanical/non-regression surfaces while YA owns the likelihood tests. |
| KD | 4.1, 4.3-4.7 | Begin 4.1 after its Phase 1 dependencies are stable. Defer power-analysis content until YA decides what ships. |
| YA + KD | 4.2 | YA leads technical examples after the Phase 2 API is frozen; KD reviews narrative, build behavior, and declared dependencies. |
| KD | 5.1-5.6 | Coordinate changed regression documentation with YA before regenerating final Rd files. |
| YA + KD | 6.1-6.8 | KD leads checks and `cran-comments.md`; YA owns maintainer-address confirmation, release tag, submission, and reviewer replies. Do not begin until YA signs off 2.10-2.15. |

Immediate parallel start:

1. YA decides the simulation scope; decisions 2.2 and 2.9 are already approved.
2. KD continues 1.1-1.17. YA starts Phase 2 only when implementation is explicitly requested.
3. KD creates the testthat scaffold (3.1) early so both contributors can add tests without colliding.

## Open questions to settle between the two of you

1. **Breaking the formula interface (2.2; approved 2026-08-25).** Move `A` and `Y`
   out of the formula RHS into explicit arguments. This intentionally breaks the old
   positional interface; update examples, tests, vignettes, and `NEWS.md` when implemented.
2. **Deprecate dotted names (2.9; approved 2026-08-25).** Introduce underscore names
   and retain the old dotted names as `.Deprecated()` wrappers for one release. Proceed
   with the real S3 methods planned in 2.5. No implementation has started.
3. **How much of `sim.*` ships in 0.1.0?** Six simulation functions with heavy runtimes
   are a CRAN check-time liability. Shipping `bc.est`, `cmBound`, `cmreg`, `cmreg.p`, and
   the predict functions in 0.1.0 and holding the power-analysis suite for 0.2.0 is a
   defensible way to get on CRAN sooner. Yuki's call.
4. **Paper citation (resolved 2026-08-25).** Atsusaka, Y., & Stevenson, R. T.
   (2023). "A Bias-Corrected Estimator for the Crosswise Model with Inattentive
   Respondents." *Political Analysis*, 31(1), 134-148.
   <https://doi.org/10.1017/pan.2021.43>.
