# cWise 0.1.0

## Breaking changes

- `cmreg()` now takes the anchor response through an explicit `anchor` argument;
  remove the anchor from the right-hand side of `formula`.
- `cmreg_p()` now takes the crosswise and anchor responses through explicit
  `crosswise` and `anchor` arguments; remove both from the right-hand side of
  `formula`.
- Underscore names are now canonical: `bc_est`, `cmreg_p`, `cmpredict_p`,
  `sim_cwdata`, `sim_estimates`, `sim_power`, and `sim_power_N`. The former
  dotted names remain as deprecated wrappers for one release.
- The canonical simulation functions now use `prevalence` rather than `pi`.
  Deprecated dotted wrappers translate their former `pi =` argument during the
  transition.

## CRAN submission preparation

- Updated package authorship metadata so that Yuki Atsusaka is the maintainer and
  Kolbe Dumas and Randy T. Stevenson are authors.
- Replaced the pre-publication citation with Atsusaka and Stevenson (2023),
  *Political Analysis*, <doi:10.1017/pan.2021.43>.
- Declared direct package dependencies and vignette-building dependencies required
  for CRAN checks.
- Revised the package title and description for CRAN metadata requirements.
- Removed legacy duplicate power-analysis code. The undocumented internal
  `sim.curve()` helper was already broken and is no longer included.
- `bc_est()` now accepts an optional `seed` for reproducible bootstrapping and
  restores the caller's random-number state afterward.
- Simulation functions now accept `verbose` to silence progress bars and
  messages. `cmBound()` examples are runnable without attaching `ggplot2`, and
  use current ggplot2 line-width syntax.
- Excluded development-only files from source builds and moved the Appendix C5
  reference PDF to `inst/doc_ref/` so it ships as a package reference.
