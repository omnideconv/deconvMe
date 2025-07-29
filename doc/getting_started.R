## ----eval=FALSE---------------------------------------------------------------
#  # install the `pak` package manager
#  install.packages("pak")
#  
#  pak::pkg_install("omnideconv/methyldeconv")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(methyldeconv)
library(minfi)
library(minfiData)

# use example data from Minfi
methyl_set <- minfiData::MsetEx
ratio_set <- minfi::ratioConvert(methyl_set)
beta_matrix <- minfi::getBeta(ratio_set)

# run EpiDISH for deconvolution of example data
result <- methyldeconv::deconvolute(methyl_set = methyl_set, method = 'epidish')

## -----------------------------------------------------------------------------
knitr::kable(head(result), caption = "Estimated cell-type fractions for the first few samples")

