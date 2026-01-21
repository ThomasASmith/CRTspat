
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CRTspat

<!-- badges: start -->

[![R-CMD-check](https://github.com/ThomasASmith/CRTspat/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ThomasASmith/CRTspat/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

[`CRTspat`](https://thomasasmith.github.io/index.html) is an R-package
to provide easily accessible R functions for design and analysis of
cluster randomized trials (CRTs), where geographically structured
spillover is anticipated and where geolocations are available. It
includes functions for analysing baseline data, for defining clusters by
algorithm, for power and sample size calculation, and for analysis of
trial outcomes. It is designed to function with outcomes that are
proportions, count data, or continuous variables.

The package was developed with CRTs of malaria interventions in mind,
where the spillover is assumed to arise as a result of mosquito
movement, with mosquito dispersal approximated with a simple diffusion
model. This does not preclude its use in other fields of research. The
anticipated use cases are described in the vignettes (articles)

The package builds on the work of [Multerer *et al.*
(2021a)](https://link.springer.com/article/10.1186/s13063-021-05543-8),
[Multerer *et al.*
(2021b)](https://link.springer.com/article/10.1186/s12936-021-03924-7)
and [Anaya-Izquierdo &
Alexander(2021)](https://onlinelibrary.wiley.com/doi/full/10.1111/biom.13316).

## Code of Conduct

Please note that the CRTspat project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Installation

`CRTspat` is on [CRAN](https://CRAN.R-project.org/package=CRTspat) and
on [github](https://github.com/ThomasASmith/CRTspat/). It can be
installed from CRAN.
