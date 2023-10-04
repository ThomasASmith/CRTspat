
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CRTspat

<!-- badges: start -->

[![R-CMD-check](https://github.com/ThomasASmith/CRTspat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ThomasASmith/CRTspat/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

[`CRTspat`](https://thomasasmith.github.io/articles/CRTspat.html) is an
R-package to provide easily accessible R functions for design and
analysis of cluster randomized trials (CRTs), where geographically
structured contamination is anticipated and where geolocations are
available. It includes functions for analysing baseline data, for
defining clusters by algorithm, for power and sample size calculation,
and for analysis of trial outcomes.

The package was developed with CRTs of malaria interventions in mind,
where the contamination is assumed to arise as a result of mosquito
movement, with mosquito dispersal approximated with a simple diffusion
model.

## Code of Conduct

Please note that the CRTspat project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Installation

`CRTspat` is on [github](https://github.com/ThomasASmith/CRTspat/) and
can be installed using:

`install.packages("devtools")` `library(devtools)`
`install_github("thomasasmith/CRTspat")`

`CRTspat` has also been submitted to CRAN.
