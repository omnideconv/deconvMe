# suppressMessages(library(methylDeconv))
# suppressMessages(library(minfiData))
# suppressMessages(library(minfi))
# 
# methyl_set <- minfiData::MsetEx
# result <- methyldeconv::deconvolute(methyl_set =  methyl_set, method = 'epidish')

suppressMessages(library(methyldeconv))

# EpiDISH signature matrix

test_that("get_epidish_signature_matrix returns a data frame with CpGs column", {
  skip_if_not_installed("EpiDISH")
  sig <- methyldeconv::get_epidish_signature_matrix("blood")
  expect_s3_class(sig, "data.frame")
  expect_true(nrow(sig) > 0)
  expect_true(ncol(sig) > 0)
  expect_true("CpGs" %in% colnames(sig))
})

# Houseman signature matrix

test_that("get_houseman_signature_matrix returns a data frame with CpGs column", {
  skip_if_not_installed("FlowSorted.Blood.EPIC")
  sig <- methyldeconv::get_houseman_signature_matrix()
  expect_s3_class(sig, "data.frame")
  expect_true(nrow(sig) > 0)
  expect_true(ncol(sig) > 0)
  expect_true("CpGs" %in% colnames(sig))
})

# MethAtlas signature matrix

test_that("get_methatlas_signature_matrix returns a data frame with expected structure", {
  ref_path <- system.file("reference_atlas.csv", package = "methyldeconv")
  sig <- methyldeconv::get_methatlas_signature_matrix(ref_path)
  expect_s3_class(sig, "data.frame")
  expect_true(nrow(sig) > 0)
  expect_true(ncol(sig) > 0)
  expect_true("CpGs" %in% colnames(sig))
})

# MethylResolver signature matrix

test_that("get_methylresolver_signature_matrix returns a data frame with CpGs column", {
  skip_if_not_installed("MethylResolver")
  sig <- methyldeconv::get_methylresolver_signature_matrix()
  expect_s3_class(sig, "data.frame")
  expect_true(nrow(sig) > 0)
  expect_true(ncol(sig) > 0)
  expect_true("CpGs" %in% colnames(sig))
})

# MethylCC signature matrix

test_that("get_methylcc_signature_matrix returns a data frame with region columns", {
  skip_if_not_installed("FlowSorted.Blood.450k")
  skip_if_not_installed("methylCC")
  sig <- methyldeconv::get_methylcc_signature_matrix()
  expect_s3_class(sig, "data.frame")
  expect_true(nrow(sig) > 0)
  expect_true(ncol(sig) > 0)
  expect_true(all(c("seqnames", "start", "end") %in% colnames(sig)))
})

