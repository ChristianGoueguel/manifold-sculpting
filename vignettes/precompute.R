# Precomputes the vignettes whose computations are too long to run when the
# package is checked: the comparison with other methods (about an hour) and
# the NIR spectra of meat (about ten minutes).
#
# Run from the root of the package source, with the packages listed in each
# vignette installed. The code is in the .Rmd.orig files; knitting them writes
# the .Rmd files, with the results and the figures (comparison-*.png and
# spectroscopy-*.png) in place, so that the vignettes build without running
# the computations. Pass the names of the vignettes to knit only some of them:
#
#   Rscript vignettes/precompute.R spectroscopy

vignettes <- commandArgs(trailingOnly = TRUE)
if (length(vignettes) == 0) vignettes <- c("comparison", "spectroscopy")

old <- setwd("vignettes")
for (v in vignettes) {
  knitr::knit(paste0(v, ".Rmd.orig"), output = paste0(v, ".Rmd"))
}
setwd(old)
