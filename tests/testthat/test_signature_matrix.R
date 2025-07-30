suppressMessages(library(deconvMe))
options(matrixStats.useNames.NA = "deprecated")

methyl_set <- minfiData::MsetEx
ratio_set <- minfi::ratioConvert(methyl_set)
beta_matrix <- minfi::getBeta(ratio_set)


test_that("Houseman, EpiDISH, MethAtlas, and MethylResolver signatures have same design",{
    houseman_sig <- deconvMe::get_houseman_signature_matrix()
    epidish_sig <- deconvMe::get_epidish_signature_matrix("blood")
    methatlas_sig <- deconvMe::get_methatlas_signature_matrix()
    methylresolver_sig <- deconvMe::get_methylresolver_signature_matrix()
    
    expect_true('CpGs' %in% colnames(houseman_sig))
    expect_true('CpGs' %in% colnames(epidish_sig))
    expect_true('CpGs' %in% colnames(methatlas_sig))
    expect_true('CpGs' %in% colnames(methylresolver_sig))
    
    expect_true(is.data.frame(houseman_sig))
    expect_true(is.data.frame(epidish_sig))
    expect_true(is.data.frame(methatlas_sig))
    expect_true(is.data.frame(methylresolver_sig))
    
    expect_true(nrow(houseman_sig) > 0)
    expect_true(nrow(epidish_sig) > 0)
    expect_true(nrow(methatlas_sig) > 0)
    expect_true(nrow(methylresolver_sig) > 0)
})

test_that("run_epidish works with Houseman signature as external reference", {
  houseman_sig <- deconvMe::get_houseman_signature_matrix()
  res <- suppressWarnings(deconvMe::run_epidish(beta_matrix, reference = houseman_sig))
  
  expect_true(is.list(res))
  expect_true("estF" %in% names(res))
  expect_true(nrow(res$estF) > 0)
  expect_true(ncol(res$estF) > 0)
})

test_that("run_methatlas works with Houseman signature as external reference", {
  houseman_sig <- deconvMe::get_houseman_signature_matrix()
  res <- suppressWarnings(deconvMe::run_methatlas(beta_matrix, reference = houseman_sig))
  
  expect_true(is.matrix(res) || is.data.frame(res))
  expect_true(nrow(res) > 0)
  expect_true(ncol(res) > 0)
})

test_that("run_methylresolver works with Houseman signature as external reference", {
  houseman_sig <- deconvMe::get_houseman_signature_matrix()
  res <- suppressWarnings(deconvMe::run_methylresolver(beta_matrix, reference = houseman_sig, alpha = 1))
  
  expect_true(is.list(res))
  expect_true("result_fractions" %in% names(res))
  expect_true(nrow(res$result_fractions) > 0)
  expect_true(ncol(res$result_fractions) > 0)
})

test_that("run_methylresolver works with Methatlas signature as external reference", {
  methatlas_sig <- deconvMe::get_methatlas_signature_matrix()
  res <- suppressWarnings(deconvMe::run_methylresolver(beta_matrix, reference = methatlas_sig, alpha = 1))
  
  expect_true(is.list(res))
  expect_true("result_fractions" %in% names(res))
  expect_true(nrow(res$result_fractions) > 0)
  expect_true(ncol(res$result_fractions) > 0)
})

