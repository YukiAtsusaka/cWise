# cWise: A (Cross)Wise Method to Analyze Sensitive Survey Questions

This package implements the methods proposed by Atsusaka and Stevenson (2020) "Bias-Corrected Crosswise Estimators for Sensitive Inquiries" (Working Paper). Our workhorse function is `bc.est` which generates a bias-corrected crosswise estimate of the proportion of individuals with sensitive attributes. `cmBound` produces applies our sensitivity analysis to crosswise data without the anchor question. `cmreg` and `cmreg.p` implement crosswise regressions in which the latent sensitive trait can be used as an outcome or as a predictor, respectively. A simulated crosswise data is saved as `cmdata`.



## Instllation
To install the latest development version of `overviewR` directly from
[GitHub](https://github.com/YukiAtsusaka/cWise) use:

``` r
library(devtools)
devtools::install_github("YukiAtsusaka/cWise")
```

## Example

First, load the package.

``` r
library(cWise)
```

The following examples use a toy data set (`cmdata`) that comes with
the package. This data contains artificially generated information in a survey using the crosswise model.

``` r
data(cmdata)
head(cmdata)
#>   ccode year month      gdp population
#> 1   RWA 1990   Jan 24180.77  14969.988
#> 2   RWA 1990   Feb 23650.53  11791.464
#> 3   RWA 1990   Mar 21860.14  30047.979
#> 4   RWA 1990   Apr 20801.06  19853.556
#> 5   RWA 1990   May 18702.84   5148.118
#> 6   RWA 1990   Jun 30272.37  48625.140
```
