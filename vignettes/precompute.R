# Precomputes the comparison vignette, which takes about an hour.
#
# Run from the root of the package source, with the packages listed in
# vignettes/comparison.Rmd.orig installed. The code is in comparison.Rmd.orig;
# knitting it writes comparison.Rmd, with the results and the figures
# (comparison-*.png) in place, so that the vignette builds without running
# the comparison.

old <- setwd("vignettes")
knitr::knit("comparison.Rmd.orig", output = "comparison.Rmd")
setwd(old)
