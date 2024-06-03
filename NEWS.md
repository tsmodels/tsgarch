# tsgarch 1.0.3

* Added the log-likelihood vector to the returned fitted and filtered object
as this will be needed in the calculation of the standard errors in the 
multivariate GARCH package.
* Added an extra option to the estimation method which adds the tmb object
to the returned estimation object. This can then be used to directly vary
the parameters and quickly extract information from the filtration. This
could be also achieved by the tsfilter method with a specification object
but is much slower. Use case is for the multivariate GARCH partitioned 
hessian calculation.

# tsgarch 1.0.2

* Moved a unit tests back to original folder and added a tolerance
to the expectation per CRAN maintainers directions.

# tsgarch 1.0.1

* Moved a couple of unit tests to other folder to avoid checking on CRAN
since the M1 mac had different rounding errors than other architectures
for simulation tests.


# tsgarch 1.0.0

* Initial CRAN submission.
* Changes to initialization of recursion in the ARCH equation to be more consistent
with the literature. This leads to a more than doubling in the accuracy against the
FCP benchmark.
* Correction to EGARCH forecast to account for log bias.
* Added a couple more data series for benchmarking.
* Added a pdf vignette.
* Added demo html vignettes.
* Added extensive unit tests.

