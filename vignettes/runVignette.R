#
# since the vignette takes a long time to run, it is not run with every build
# it has been following the suggestion here:
# as suggested here https://ropensci.org/blog/2019/12/08/precompute-vignettes/
# To run the vignettes, execute the code from the vignettes via:


library(CRTspat)
knitr::opts_chunk$set(error=FALSE)
knitr::knit("vignettes/CRTspat.Rmd.orig", output = "vignettes/CRTspat.Rmd")
knitr::knit("vignettes/Usecase1.Rmd.orig", output = "vignettes/Usecase1.Rmd")
knitr::knit("vignettes/Usecase2.Rmd.orig", output = "vignettes/Usecase2.Rmd")
knitr::knit("vignettes/Usecase3.Rmd.orig", output = "vignettes/Usecase3.Rmd")
rmarkdown.html_vignette.check_title = FALSE
devtools::install(build_vignettes = TRUE)


# to build package website
usethis::use_pkgdown()
pkgdown::build_site()

