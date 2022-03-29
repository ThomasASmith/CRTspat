
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CRTspillover

<!-- badges: start -->

<!-- badges: end -->

The goal of CRTspillover is to provide easily accessible R functions for
designing cluster randomized trials (CRTs), where geographically
structured contamination is anticipated. The package was developed using
examples from CRTs of malaria interventions, where the contamination is
assumed to arise as a result of mosquito movement, with mosquito
dispersal approximated with a simple diffusion model.

## Installation

You can install the released version of CRTspillover from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("CRTspillover")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(CRTspillover)
## Example generation of default assignments for a test dataset
testDesign = Design_CRT() 
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

    {r cars}
    summary(cars)

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
