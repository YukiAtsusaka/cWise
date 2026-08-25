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
