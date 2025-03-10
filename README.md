# methyldeconv

<!-- badges: start -->

[![R-CMD-check](https://github.com/omnideconv/methyldeconv/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/omnideconv/methyldeconv/actions/workflows/R-CMD-check.yml) [![Codecov test coverage](https://codecov.io/gh/omnideconv/methyldeconv/branch/main/graph/badge.svg)](https://app.codecov.io/gh/omnideconv/methyldeconv?branch=main)

<!-- badges: end -->

Ever wanted to apply cell-type deconvolution on your DNA methylation data but could not decide which method to use? Here is **methyldeconv** to help your needs!

This package integrates unified access to five reference-based cell-type deconvolution methods that can directly be applied to Illumina array data (450k, EPIC arrays) or bisulfite sequencing data (RRBS, WGBS).

The included methods are:

| method  | license | citation|
|------------------------|------------------------|------------------------|
| [EpiDISH](https://bioconductor.org/packages/release/bioc/html/EpiDISH.html)                                                   | GPL-2   | Teschendorff, A.E., Breeze, C.E., Zheng, S.C. *et al.* A comparison of reference-based algorithms for correcting cell-type heterogeneity in Epigenome-Wide Association Studies. *BMC Bioinformatics* **18**, 105 (2017). <https://doi.org/10.1186/s12859-017-1511-5>|
| [Houseman (Flow.Sorted.Blood)](https://www.bioconductor.org/packages/release/data/experiment/html/FlowSorted.Blood.EPIC.html) | GPL-3   | Houseman, E.A., Accomando, W.P., Koestler, D.C. *et al.* DNA methylation arrays as surrogate measures of cell mixture distribution. *BMC Bioinformatics* **13**, 86 (2012). <https://doi.org/10.1186/1471-2105-13-86> Koestler, D.C., Jones, M.J., Usset, J. et al. Improving cell mixture deconvolution by identifying optimal DNA methylation libraries (IDOL).BMC Bioinformatics 17, 120 (2016). <https://doi.org/10.1186/s12859-016-0943-7> Salas, L.A., Koestler, D.C., Butler, R.A. et al. An optimized library for reference-based deconvolution of whole-blood biospecimens assayed using the Illumina HumanMethylationEPIC BeadArray. Genome Biol 19, 64 (2018). <https://doi.org/10.1186/s13059-018-1448-7> |
|[MethAtlas](https://github.com/nloyfer/meth_atlas)|[Yissum Software Research License](https://github.com/nloyfer/meth_atlas/blob/master/LICENSE.md) "This software is distributed under the **Yissum Software Research License**, which permits use, modification, and distribution for **research purposes only**. Commercial use requires a separate license from Yissum (software\@yissum.co.il)."|Moss, J., Magenheim, J., Neiman, D. *et al.* Comprehensive human cell-type methylation atlas reveals origins of circulating cell-free DNA in health and disease. *Nat Commun* **9**, 5068 (2018). <https://doi.org/10.1038/s41467-018-07466-6>|
|[methylCC](https://github.com/stephaniehicks/methylCC)|GPL-3|Hicks, S.C., Irizarry, R.A. methylCC: technology-independent estimation of cell type composition using differentially methylated regions. *Genome Biol* **20**, 261 (2019). <https://doi.org/10.1186/s13059-019-1827-8>|
|[methylResolver](https://github.com/darneson/MethylResolver)|GPL-3|Arneson, D., Yang, X. & Wang, K. MethylResolver—a method for deconvoluting bulk DNA methylation profiles into known and unknown cell contents. *Commun Biol* **3**, 422 (2020). <https://doi.org/10.1038/s42003-020-01146-2>|

## Installation

You can install methyldeconv from [GitHub](https://github.com/), we recommend to use the [pak](https://github.com/r-lib/pak) package manager:

``` r
# install the `pak` package manager
install.packages("pak")

pak::pkg_install("omnideconv/methyldeconv")
```

## Example

methyldeconv can either be applied directly to a methylSet from the minfi package, or you can apply each method separately on a beta matrix with Illumina CpG IDs.

Both cases will be demonstrated here using example data from minfi:

``` r
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
```

With methyldeconv you can either get the original result object of each respective method (`methyldeconv::run_XXX`) or a table with cell-type fractions that has a unified format between methods (`methyldeconv::deconvolute`).The unified results can be visualized using the methyldeconv functions `results_barplot()` or `results_boxplot()`.

``` r
methyldeconv::results_barplot(result)

methyldeconv::results_boxplot(result)
```

Results from a run with more than one method (`methyldeconv::deconvolute_combined`) can be easily visualized as well with `results_aggregated_boxplot`:

```r
methyldeconv::results_aggregated_boxplot(result_multiple)
```


## Dependencies

A full list of dependencies can be displayed with `pak`:
``` r
pak::pkg_deps_tree("omnideconv/methyldeconv")
```