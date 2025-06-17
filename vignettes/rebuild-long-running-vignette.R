old_wd <- getwd()

setwd("vignettes/")
knitr::knit("blindrecalc_paper.Rmd.orig", output = "blindrecalc_paper.Rmd")
knitr::purl("blindrecalc_paper.Rmd.orig", output = "blindrecalc_paper.R")

setwd(old_wd)
