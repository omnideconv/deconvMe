suppressMessages(library(deconvMe))
suppressMessages(library(minfiData))
suppressMessages(library(minfi))

methyl_set <- minfiData::MsetEx
ratio_set <- minfi::ratioConvert(methyl_set)
beta_matrix <- minfi::getBeta(ratio_set)



test_that("EpiDISH works", {
  epidish_res <- deconvMe::run_epidish(beta_matrix = beta_matrix, mode='RPC')$estF
  check_result <- as.matrix(read.csv("test_results/epidish.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = epidish_res,
    expected = check_result, tolerance = 1e-3
  )
})


test_that("methylCC works", {
  methylcc_res <- as.matrix(deconvMe::run_methylcc(methyl_set = methyl_set))
  check_result <- as.matrix(read.csv("test_results/methylcc.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = methylcc_res,
    expected = check_result, tolerance = 1e-3
  )
})

test_that("Houseman works", {
  flowSorted_res <- deconvMe::run_houseman(methyl_set = methyl_set)$prop
  check_result <- as.matrix(read.csv("test_results/houseman.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = flowSorted_res,
    expected = check_result, tolerance = 1e-3
  )
})

test_that("MethylResolver works", {
  methylResolver_res <- as.matrix(deconvMe::run_methylresolver(beta_matrix = beta_matrix, alpha = 1)$result_fractions)
  check_result <- as.matrix(read.csv("test_results/methylresolver.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = methylResolver_res,
    expected = check_result, tolerance = 1e-3
  )
})

test_that("MethAtlas works", {
  meth_atlas_res <- deconvMe::run_methatlas(beta_matrix = beta_matrix)
  check_result <- as.matrix(read.csv("test_results/methatlas.csv",
                                     row.names = 1,
                                     check.names = FALSE
  ))
  expect_equal(
    info = "deconvolution result is correct", object = meth_atlas_res,
    expected = check_result, tolerance = 1e-3
  )
})

test_that("Main function deconvolute works (tested on epidish)", {
  res <- deconvMe::deconvolute(methyl_set = methyl_set, method = 'epidish')
  res <- res |> dplyr::select(order(colnames(res)))
  check_result <- read.csv("test_results/epidish.csv",
                                     row.names = 1,
                                     check.names = FALSE
  )
  check_result <- check_result |> dplyr::select(order(colnames(check_result)))
  expect_equal(
    info = "deconvolution result is correct", object = res,
    expected = check_result, tolerance = 1e-3
  )
})

test_that("Running multiple methods works", {
  res <- deconvMe::deconvolute_combined(methyl_set = methyl_set,
                                            methods = c('epidish','houseman'),
                                            array = '450k')

  methods <- unique(res$method)
  expect_equal(
    object = methods, expected = c('epidish','houseman','aggregated'),
    info = 'Correct method outputs and aggregated outputs are present.'
  )
})

