## Resubmission: 2026-09-17

This resubmission addresses Leonore Hochhauser's review of the initial cWise
0.1.0 submission:

* The DESCRIPTION reference uses `Atsusaka and Stevenson (2023)
  <doi:10.1017/pan.2021.43>`. This is already present in the current repository
  source and has been verified in the newly built submission archive.
* Replaced `dontrun` with `donttest` in the four lengthy simulation examples:
  `sim_cwdata`, `sim_estimates`, `sim_power`, and `sim_power_N`, including their
  generated Rd files. Timed executions took approximately 7, 28, 158, and
  86 seconds, respectively, on the local macOS host. Short examples remain
  unwrapped.
* No estimator, API, authorship, or vignette-source changes were made for this
  review response. The version remains 0.1.0 because it has not yet been
  accepted on CRAN.

## Current resubmission validation

* macOS 26.5.2, Apple Silicon, R 4.4.2 (2024-10-31):
  `devtools::check(args = c("--as-cran", "--run-donttest"))`, which ran
  `R CMD check --as-cran --run-donttest --no-manual`.
* Results: **0 errors, 0 warnings, and 1 note**. The host could not verify the
  current time online.
* All examples passed, including the four `donttest` blocks. All 52 test
  assertions passed with no test warnings, failures, or skips. The current
  combined vignette built and rebuilt successfully and is included in the archive.
* CRAN incoming and remote incoming checks were disabled by the local devtools
  check environment. This result does not establish those gates or a fresh
  Windows/Linux check result.
* A separate default `R CMD Rd2pdf` attempt was blocked because the local TeX
  installation lacks `inconsolata.sty`. Rd syntax, metadata, cross-references,
  usage, and content checks passed. A default PDF manual check remains pending
  on a host with that font package installed.

## Previous submission checks (2026-09-08)

The following records describe the previous submission, not a fresh check of
this resubmission.

### Test environments

* Windows 11 x64 (build 26200), R 4.4.3 (2025-02-28):
  `R CMD check --as-cran --no-manual`.
* Windows r-devel: submitted to win-builder.r-project.org on 2026-09-08;
  report pending at `atsusaka@uh.edu`.
* Linux and macOS: R-hub check pending. The official GitHub Actions workflow
  (`.github/workflows/rhub.yaml`) is prepared locally and must be committed and
  pushed before R-hub can dispatch the checks.

### R CMD check results

The local check completed with **0 errors, 0 warnings, and 4 notes**:

* The CRAN incoming and URL checks could not run because this Windows host's
  TLS credentials are unavailable (`schannel` error 35). The DOI, GitHub, and
  package website URLs are valid and must be rechecked from a host with working
  HTTPS before submission.
* The host could not verify the current time online.
* Pandoc was not discoverable through the standard R check environment; the
  package vignettes were nevertheless built and rebuilt successfully using the
  installed Quarto Pandoc.
* `cmreg_p` and `cmpredict_p` examples each took slightly more than five
  seconds. They remain well within CRAN's total check-time budget.

This is a new submission. Update this file with the win-builder and R-hub
results before running `devtools::release()`.
