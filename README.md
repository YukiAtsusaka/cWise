# cWise: A (Cross)Wise Method to Analyze Sensitive Survey Questions 

<!-- badges: start -->

[![R
badge](https://img.shields.io/badge/Build%20with-🔥%20and%20R-blue)](https://github.com/cosimameyer/overviewR)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/cWise)](https://cran.r-project.org/package=cWise)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)　<img src='man/figures/lisafotios.jpg' align="right" height="200" />
<!-- [![Rdoc](https://www.rdocumentation.org/badges/version/overviewR)](https://www.rdocumentation.org/packages/overviewR) -->
<!-- [![metacran downloads](https://cranlogs.r-pkg.org/badges/overviewR)](https://cran.r-project.org/package=overviewR) -->
<!-- [![cran checks](https://cranchecks.info/badges/summary/overviewR)](https://cran.r-project.org/web/checks/check_results_overviewR.html) -->
<!-- [![](https://cranlogs.r-pkg.org/badges/version/overviewR)](https://www.r-pkg.org/badges/version/overviewR) -->
<!-- [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) -->
<!-- [![Last-changedate](https://img.shields.io/badge/last%20change-2020--07--13-green.svg)](/commits/master) -->
<!-- badges: end -->



This R package implements the methods proposed by Atsusaka and Stevenson (2020) ["Bias-Corrected Crosswise Estimators for Sensitive Inquiries"](https://github.com/YukiAtsusaka/working-paper/blob/master/WP_BiasCorrectedCM.pdf). Our workhorse function is `bc.est` which generates a bias-corrected crosswise estimate of the proportion of individuals with sensitive attributes. `cmBound` applies our sensitivity analysis to crosswise data without the anchor question. `cmreg` and `cmreg.p` implement crosswise regressions in which the latent sensitive trait can be used as an outcome or as a predictor, respectively. A simulated crosswise data is saved as `cmdata`.

<details>
<summary>Cite this software✒️</summary>

@Manual{,
    title = {cWise: A (Cross)Wise Method to Analyze Sensitive Survey Questions},
    author = {Yuki Atsusaka},
    year = {2020},
    note = {R package version 0.0.0},
    url = {https://CRAN.R-project.org/package=cWise},
  }
</details>

## Instllation
To install the latest development version of `cWise` directly from
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

#>   Y A  p      p.prime
#> 1 1 1 0.15    0.15
#> 2 0 0 0.15    0.15
#> 3 0 0 0.15    0.15
#> 4 1 0 0.15    0.15
#> 5 0 1 0.15    0.15
#> 6 1 1 0.15    0.15
```

Here, `Y` is a binary response in the crosswise question (TRUE-TRUE or FALSE-FALSE), and `A` is a binary response in the anchor question. 
`p` and `p.prime` are auxiliary probabilities in the crosswise and anchor questions, respectively.
While researchers can simply remember the values of `p` and `p.prime` (without including them in data), we include them for an illustrative purpose.
Note that `P` and `p.prime` must be (by construction) identical for all respondents.


## `bc.est`: Estimate the Prevalence of Sensitive Attributes
Generate a bias-corrected crosswise estimate using a crosswise data.

```r
prev <- bc.est(Y=Y, A=A, p=0.15, p.prime=0.15, data=cmdata)
prev

#> $Results
#>                 Estimate Std. Error 95%CI(Low) 95%CI(Up)
#> Naive Crosswise   0.2331     0.0126     0.2083    0.2579
#> Bias-Corrected    0.1547     0.0182     0.1193    0.1845
#> 
#> $Stats
#>  Attentive Rate Sample Size
#>          0.7729        2000
```


Sample weights can be easily incorporated in our bias-corrected estimator.

```r
bc.est(Y=cross, A=anchor, p=p.cross, p.prime=p.anchor, w=sweight, data=cmdata)
```



##  `cmBound`: Apply a Sensitivity Analysis
Apply sensitivity analysis and generate sensitivity bounds for naive crosswise estimates.

```r
p <- cmBound(lambda.hat=0.6385, p=0.25, N=310, dq=0.073, N.dq=310)
p
```

<img src="man/figures/bounds.png" width="50%" style="display: block; margin: auto;" />

Since the output is a ggplot object, one can easily add additional information by using "+". 
For example, to add a title with a specific font:

```r
p <- p + ggtitle("Sensitivity Analysis") + 
         theme(plot.title = element_text(size=20, face="bold"))       
p         
```

<img src="man/figures/bounds2.png" width="50%" style="display: block; margin: auto;" />



##  `cmreg`: Run a Regression with the Latent Sensitive Trait as an Outcome

For an illustration, let's load and see `cmdata2` that contains the main and anchor response variables along with two covariates.

```r
data(cmdata2)
head(cmdata2)

#>   Y A female age   p   p2
#> 1 1 1      0  23 0.1 0.15
#> 2 1 1      1  31 0.1 0.15
#> 3 0 1      1  32 0.1 0.15
#> 4 1 0      1  19 0.1 0.15
#> 5 0 1      1  19 0.1 0.15
#> 6 0 1      1  25 0.1 0.15
```

To run a crosswise regression, one can specify the model by writing a formula: `Crosswise Response ~ var1 + ... + varN + Anchor Response` as follow:

```r
m <-  cmreg(Y~female+age+A, p=0.1, p2=0.15, data=cmdata2)
m

#> $Call
#> Y ~ female + age + A
#>
#> $Coefficients
#>                Estimate Std. Error
#> (intercept) -1.65085102 0.42681118
#> female       0.28162406 0.14267305
#> age          0.03264242 0.01332295
#>
#> $AuxiliaryCoef
#>               Estimate Std. Error
#> (intercept)  0.13868287 1.13470799
#> female      -0.20436648 0.41193379
#> age          0.05945707 0.03935859
```

##  `cmreg.p`: Run a Regression with the Latent Sensitive Trait as a Predictor



```r
data(cmdata3)
head(cmdata3)

#>             V Y female age A   p   p2
#> 1 -0.38350925 1      0  23 1 0.1 0.15
#> 2 -0.05965305 1      1  31 1 0.1 0.15
#> 3  0.72655660 0      1  32 1 0.1 0.15
#> 4  0.79845870 1      1  19 0 0.1 0.15
#> 5 -0.19410532 0      1  19 1 0.1 0.15
#> 6 -0.34926673 0      1  25 1 0.1 0.15
```

```r
m2 <- cmreg.p(V~age+female+Y+A, p=0.1, p2=0.15, data=cmdata3)
m2

#> $Call
#> V ~ age + female + Y + A
#>
#> $Coefficients
#>                Estimate  Std. Error
#> (intercept) 0.023517624 0.147752825
#> age         0.009642705 0.004806593
#> female      0.247318563 0.052038635
#> Y           0.985775066 0.075624764
#>
#> $AuxiliaryCoef
#>               Estimate Std. Error
#> (intercept) -1.7338112 0.40093115
#> age          0.0351617 0.01258319
#> female       0.2878427 0.13560686
#> 
#> $AuxiliaryCoef2
#>                Estimate Std. Error
#> (intercept)  0.24805722 1.06797382
#> age          0.05479166 0.03699454
#> female      -0.10749144 0.37791738
```
