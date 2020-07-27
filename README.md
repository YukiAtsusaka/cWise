# cWise: A (Cross)Wise Method to Analyze Sensitive Survey Questions

This package implements the methods proposed by Atsusaka and Stevenson (2020) "Bias-Corrected Crosswise Estimators for Sensitive Inquiries" (Working Paper). Our workhorse function is `bc.est` which generates a bias-corrected crosswise estimate of the proportion of individuals with sensitive attributes. `cmBound` applies our sensitivity analysis to crosswise data without the anchor question. `cmreg` and `cmreg.p` implement crosswise regressions in which the latent sensitive trait can be used as an outcome or as a predictor, respectively. A simulated crosswise data is saved as `cmdata`.



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
#> cross anchor    p.cross p.anchor
#> 1 1 1 0.15    0.15
#> 2 0 0 0.15    0.15
#> 3 0 0 0.15    0.15
#> 4 1 0 0.15    0.15
#> 5 0 1 0.15    0.15
#> 6 1 1 0.15    0.15
```

## `bc.est`
Generate a bias-corrected crosswise estimate using a crosswise data.

```r
bc.est(Y=cross, A=anchor, p=p.cross, p.prime=p.anchor, data=cmdata)

```

## `cmBound`
Apply sensitivity analysis and generate sensitivity bounds for naive crosswise estimates.

```r
sensitivity <- cmBound(lambda.hat=0.6385, p=0.25, N=310, dq=0.073)
```





