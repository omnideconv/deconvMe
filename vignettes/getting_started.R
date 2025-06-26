## ----eval=FALSE---------------------------------------------------------------
# # install the `pak` package manager
# install.packages("pak")
# 
# pak::pkg_install("omnideconv/methyldeconv")

## -----------------------------------------------------------------------------
library(methyldeconv)
library(minfi)
library(minfiData)

# use example data from Minfi
methyl_set <- minfiData::MsetEx
ratio_set <- minfi::ratioConvert(methyl_set)
beta_matrix <- minfi::getBeta(ratio_set)

# run EpiDISH for deconvolution of example data
result <- methyldeconv::deconvolute(methyl_set = methyl_set, method = 'epidish')

result_raw <- methyldeconv::run_epidish(beta_matrix = beta_matrix, mode='RPC')

# you can also run multiple methods at the same time and get their results + aggregated results:
result_multiple <- methyldeconv::deconvolute_combined(methyl_set = methyl_set,
                                                      methods = c('epidish','houseman'),
                                                      array = '450k')

## -----------------------------------------------------------------------------
methyldeconv::results_barplot(result)

methyldeconv::results_boxplot(result)

## -----------------------------------------------------------------------------
methyldeconv::results_aggregated_boxplot(result_multiple)

