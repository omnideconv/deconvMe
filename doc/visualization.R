## ----message=FALSE, warning=FALSE---------------------------------------------
library(methyldeconv)
library(minfi)
library(minfiData)

# Example data
methyl_set <- minfiData::MsetEx
ratio_set <- minfi::ratioConvert(methyl_set)
beta_matrix <- minfi::getBeta(ratio_set)

# Run deconvolution
result <- methyldeconv::deconvolute(methyl_set = methyl_set, method = 'epidish')
result_multiple <- methyldeconv::deconvolute_combined(methyl_set = methyl_set,
                                                      methods = c('epidish','methylcc'),
                                                      array = '450k')

## ----fig.height=4, fig.width=8------------------------------------------------
methyldeconv::results_barplot(result)

## ----fig.height=4, fig.width=8------------------------------------------------
methyldeconv::results_boxplot(result)

## ----fig.height=4, fig.width=8------------------------------------------------
methyldeconv::results_aggregated_boxplot(result_multiple)

