
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--04--21-yellowgreen.svg)](/commits/master)
[![packageversion](https://img.shields.io/badge/Package%20version-1.0.0-orange.svg?style=flat-square)](commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/tsgarch)](https://cran.r-project.org/package=tsgarch)

# tsgarch

[A new version for CRAN submission is under development with many fixes and more accurate benchmarks]

The **tsgarch** package is a partial re-implementation of
[rugarch](https://CRAN.R-project.org/package=rugarch), by the same
author, with key differences summarized below:

- it does not (yet) implement all GARCH models in rugarch. FIGARCH,
  Multiplicative Component GARCH and Realized GARCH are not currently
  implemented.
- it does not implement joint ARFIMA-GARCH estimation. The conditional
  mean equation only allows for a constant. With so many options for
  modelling the conditional mean, many of which are available in the
  **tsmodels** framework, it was decided to keep this package simpler
  and avoid the joint estimation of both conditional mean and variance
  dynamics. While the 2 step estimation approach, whereby the residuals
  of the conditional mean are passed to the variance dynamics
  estimation, may be less efficient for small sized datasets, it is
  expected to be more flexible in what can be achieved. Additionally,
  the ARCH-in-mean model is no longer available, as it was found to have
  very limited value within the **tsmodels** framework or in this
  author’s experience. A separate
  [tsarma](https://github.com/tsmodels/tsarma) package for ARMA(p,q)-X
  models is however available.
- it makes use of automatic differentiation (autodiff) during
  estimation, via the [TMB](https://CRAN.R-project.org/package=TMB)
  package. This is in line with similar approaches in other models
  written in the **tsmodels** framework. Using autodiff allows for more
  confident estimation and more accurate standard errors.
- it fully implements and correctly documents a number of sandwich
  estimators making use of the
  [sandwich](https://CRAN.R-project.org/package=sandwich) package
  framework (with methods for `bread` and `estfun` and `meat`/`meat_HAC`
  functions).
- it makes use of S3 methods and classes, abandoning the S4 approach of
  **rugarch**. Additionally, while making use of standard methods from
  the stats package, some of the methods are based on those exported
  from [tsmethods](https://CRAN.R-project.org/package=tsmethods),
  consistent with other packages in the **tsmodels** framework.

## Installation

The package can be installed from CRAN or the tsmodels github repo:

``` r
remotes::install_github("tsmodels/tsgarch", dependencies = TRUE)
```

Note, that in order to make use of symbolic output, flextable requires
[equatags](https://CRAN.R-project.org/package=equatags) to be installed
which has a dependency on
[xlst](https://CRAN.R-project.org/package=xslt) which in turn has
SystemRequirements libxslt. Therefore, if you are seeing `NA` printed in
place of symbols, then it is likely that xlst is not installed.

## Performance

Due to the heavy use of the `data.table` package for parameter tracking
and indexing, I have seen significant speed improvements by setting the
following :

``` r
options(datatable.optimize = 2L)
```
