#
# since the vignette takes a long time to run, it is not run with every build
# it has been following the suggestion here:
# as suggested here https://ropensci.org/blog/2019/12/08/precompute-vignettes/
# To run the vignettes, execute the code from the vignettes via:

knitr::knit("vignettes/CRTspillover.Rmd.orig", output = "vignettes/CRTspillover.Rmd")
knitr::knit("vignettes/Usecase1.Rmd.orig", output = "vignettes/Usecase1.Rmd")
knitr::knit("vignettes/Usecase2.Rmd.orig", output = "vignettes/Usecase2.Rmd")
knitr::knit("vignettes/Usecase3.Rmd.orig", output = "vignettes/Usecase3.Rmd")
rmarkdown.html_vignette.check_title = FALSE
devtools::install(build_vignettes = TRUE)

