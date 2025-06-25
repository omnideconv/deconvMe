library(methyldeconv)

skip_if_not_installed("minfiData")
library(minfiData)
library(minfi)

set.seed(42)
methyl_set <- minfiData::MsetEx
# Use all 6 samples and 50,000 random CpGs for a more representative subset
sample_idx <- 1:6
cpg_idx <- sample(seq_len(nrow(methyl_set)), 50000)
methyl_set_small <- methyl_set[cpg_idx, sample_idx]
ratio_set <- minfi::ratioConvert(methyl_set_small)
beta_matrix <- minfi::getBeta(ratio_set)

# Parameter variation: run_epidish

test_that("run_epidish works for all modes", {
  modes <- c("RPC", "CBS", "CP")
  for (mode in modes) {
    res <- run_epidish(beta_matrix = beta_matrix, mode = mode)
    expect_true(is.list(res))
    expect_true("estF" %in% names(res))
  }
})

# Parameter variation: run_houseman

test_that("run_houseman works for both arrays", {
  arrays <- c("450k", "EPIC")
  for (array in arrays) {
    res <- run_houseman(methyl_set = methyl_set_small, array = array)
    expect_true(is.list(res) || is.data.frame(res))
  }
})

# Parameter variation: run_methylcc

test_that("run_methylcc works for both arrays", {
  arrays <- c("450k", "EPIC")
  for (array in arrays) {
    res <- run_methylcc(methyl_set = methyl_set_small, array = array)
    expect_true(is.matrix(res) || is.data.frame(res))
  }
})

# Output structure: run_methylresolver

test_that("run_methylresolver output structure is correct", {
  res <- run_methylresolver(beta_matrix = beta_matrix, alpha = 1)
  expect_true(is.list(res))
  expect_true(all(c("result_metrics", "result_fractions", "result_absolute", "result_purity") %in% names(res)))
})

# Output structure: run_methatlas

test_that("run_methatlas output is a matrix/data.frame", {
  res <- run_methatlas(beta_matrix = beta_matrix)
  expect_true(is.matrix(res) || is.data.frame(res))
})

# Output structure: deconvolute

test_that("deconvolute output is a data.frame with expected columns", {
  res <- deconvolute(methyl_set = methyl_set_small, method = "epidish")
  colnames(res) <- methyldeconv::rename_cell_types(colnames(res))
  expect_true(is.data.frame(res))
  expect_true(all(c("B cell", "Monocyte", "Neutrophil", "NK cell", "T cell CD4+", "T cell CD8+") %in% colnames(res)))
})

# Output structure: deconvolute_combined

test_that("deconvolute_combined output contains method and aggregated", {
  res <- deconvolute_combined(methyl_set = methyl_set_small, methods = c("epidish", "houseman"), array = "450k")
  expect_true(is.data.frame(res))
  expect_true(all(c("sample", "method", "celltype", "value") %in% colnames(res)))
  expect_true("aggregated" %in% res$method)
}) 