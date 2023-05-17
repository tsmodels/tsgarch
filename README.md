
[![R-CMD-check](https://github.com/tsmodels/tsgarch/workflows/R-CMD-check/badge.svg)](https://github.com/tsmodels/tsgarch/actions)
[![Last-changedate](https://img.shields.io/badge/last%20change-2023--05--17-yellowgreen.svg)](/commits/master)
[![packageversion](https://img.shields.io/badge/Package%20version-0.1.0-orange.svg?style=flat-square)](commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/tsgarch)](https://cran.r-project.org/package=tsgarch)

# tsgarch

The **tsgarch** package is a partial re-implementation of
[rugarch](https://cran.r-project.org/web/packages/rugarch/index.html),
by the same author, with key differences summarized below:

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
  estimation, via the
  [TMB](https://cran.r-project.org/web/packages/TMB/index.html) package.
  This is in line with similar approaches in other models written in the
  **tsmodels** framework. Using autodiff allows for more confident
  estimation and more accurate standard errors.

- it fully implements and correctly documents a number of sandwich
  estimators making use of the
  [sandwich](https://cran.r-project.org/web/packages/sandwich/index.html)
  package framework (with methods for `bread` and `estfun` and
  `meat`/`meat_HAC` functions).

- it makes use of S3 methods and classes, abandoning the S4 approach of
  **rugarch**. Additionally, while making use of standard methods from
  the stats package, some of the methods are based on those exported
  from [tsmethods](https://github.com/tsmodels/tsmethods), consistent
  with other packages in the **tsmodels** framework.

The package vignette is available in the new tsmodels site
[here](https://www.nopredict.com/packages/tsgarch.html).

## Installation

The package can be installed from the tsmodels github repo. Note that
installation may take some time due to the compilation of the TMB code.

``` r
remotes::install_github("tsmodels/tsgarch", dependencies = TRUE)
```

Note, that in order to make use of symbolic output, flextable requires
[equatags](https://cran.r-project.org/web/packages/equatags/index.html)
to be installed which has a dependency on
[xlst](https://cran.r-project.org/web/packages/xslt/index.html) which in
turn has SystemRequirements libxslt. Therefore, if you are seeing `NA`
printed in place of symbols, then it is likely that xlst is not
installed.
