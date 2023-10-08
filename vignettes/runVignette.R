# To regenerate the figures change fig.keep = 'none' to fig.path=('vignettes/')
# (there is a long thread about how pkgdown doesn't seem to find figures where
# it says it does, hence the workaround with html code)

# since the vignettes take a long time to run, they are not run with every build
# as suggested here https://ropensci.org/blog/2019/12/08/precompute-vignettes/
# to run the vignettes, execute the code from the vignettes via:

library(CRTspat)
knitr::opts_chunk$set(error=FALSE)
knitr::knit("vignettes/CRTspat.Rmd.orig", output = "vignettes/CRTspat.Rmd")
knitr::knit("vignettes/Usecase1.Rmd.orig", output = "vignettes/Usecase1.Rmd")
knitr::knit("vignettes/Usecase2.Rmd.orig", output = "vignettes/Usecase2.Rmd")
knitr::knit("vignettes/Usecase3.Rmd.orig", output = "vignettes/Usecase3.Rmd")
knitr::knit("vignettes/Usecase4.Rmd.orig", output = "vignettes/Usecase4.Rmd")
knitr::knit("vignettes/Usecase5.Rmd.orig", output = "vignettes/Usecase5.Rmd")
knitr::knit("vignettes/Usecase6.Rmd.orig", output = "vignettes/Usecase6.Rmd")
knitr::knit("vignettes/Usecase7.Rmd.orig", output = "vignettes/Usecase7.Rmd")
knitr::knit("vignettes/Usecase8.Rmd.orig", output = "vignettes/Usecase8.Rmd")
knitr::knit("vignettes/Usecase9.Rmd.orig", output = "vignettes/Usecase9.Rmd")
#knitr::knit("vignettes/Usecase11.Rmd.orig", output = "vignettes/Usecase11.Rmd")

rmarkdown.html_vignette.check_title = FALSE
devtools::install(build_vignettes = TRUE)

detach("package:CRTspat", unload = TRUE)
# to build package website
usethis::use_pkgdown()
pkgdown::build_site()

# To write pdf manual
shell('R CMD Rd2pdf . --output=man/figures/manual.pdf --force --no-preview')
