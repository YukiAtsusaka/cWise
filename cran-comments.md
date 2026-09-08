## Test environments

* Windows 11 x64 (build 26200), R 4.4.3 (2025-02-28):
  `R CMD check --as-cran --no-manual`.
* Windows r-devel: submitted to win-builder.r-project.org on 2026-09-08;
  report pending at `atsusaka@uh.edu`.
* Linux and macOS: R-hub check pending. The official GitHub Actions workflow
  (`.github/workflows/rhub.yaml`) is prepared locally and must be committed and
  pushed before R-hub can dispatch the checks.

## R CMD check results

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
