# suppressMessages(library(methylDeconv))
# suppressMessages(library(minfiData))
# suppressMessages(library(minfi))
# 
# methyl_set <- minfiData::MsetEx
# result <- methyldeconv::deconvolute(methyl_set =  methyl_set, method = 'epidish')

suppressMessages(library(methyldeconv))

# ----------------------
# EpiDISH tests
# ----------------------

test_that("get_epidish_signature_matrix returns a data frame with CpGs column", {
  skip_if_not_installed("EpiDISH")
  sig <- methyldeconv::get_epidish_signature_matrix("blood")
  expect_s3_class(sig, "data.frame")
  expect_true(nrow(sig) > 0)
  expect_true(ncol(sig) > 0)
  expect_true("CpGs" %in% colnames(sig))
})

test_that("run_epidish supports cpg_subset argument", {
  skip_if_not_installed("EpiDISH")
  sig <- methyldeconv::get_epidish_signature_matrix("blood")
  subset_cpgs <- sig$CpGs[1:10]
  beta <- matrix(runif(10 * 3), nrow = 10, ncol = 3)
  rownames(beta) <- subset_cpgs
  colnames(beta) <- paste0("Sample", 1:3)
  res <- methyldeconv::run_epidish(beta, reference = "blood", cpg_subset = subset_cpgs)
  expect_true(all(rownames(res$dataREF) %in% subset_cpgs))
})

# ----------------------
# CpG overlap warning/error tests
# ----------------------

test_that("run_epidish warns and proceeds if some CpGs are missing, errors if all are missing", {
  skip_if_not_installed("EpiDISH")
  sig <- methyldeconv::get_epidish_signature_matrix("blood")
  subset_cpgs <- sig$CpGs[1:5]
  missing_cpgs <- paste0("cg_missing_", 1:5)
  beta <- matrix(runif(10 * 3), nrow = 10, ncol = 3)
  rownames(beta) <- c(subset_cpgs, missing_cpgs)
  colnames(beta) <- paste0("Sample", 1:3)
  expect_warning(
    res <- methyldeconv::run_epidish(beta, reference = "blood", cpg_subset = c(subset_cpgs, missing_cpgs)),
    "not present in the reference matrix"
  )
  expect_error(
    methyldeconv::run_epidish(beta, reference = "blood", cpg_subset = missing_cpgs),
    "None of the specified cpg_subset CpGs are present"
  )
})

# ----------------------
# Houseman tests
# ----------------------

test_that("get_houseman_signature_matrix returns a data frame with CpGs column", {
  skip_if_not_installed("FlowSorted.Blood.EPIC")
  sig <- methyldeconv::get_houseman_signature_matrix()
  expect_s3_class(sig, "data.frame")
  expect_true(nrow(sig) > 0)
  expect_true(ncol(sig) > 0)
  expect_true("CpGs" %in% colnames(sig))
})

test_that("run_houseman supports cpg_subset argument", {
  skip_if_not_installed("FlowSorted.Blood.EPIC")
  skip_if_not_installed("minfiData")
  skip_if_not_installed("minfi")
  sig <- methyldeconv::get_houseman_signature_matrix()
  subset_cpgs <- sig$CpGs[1:10]
  methyl_set <- minfiData::MsetEx
  res <- methyldeconv::run_houseman(methyl_set, cpg_subset = subset_cpgs)
  expect_true(is.data.frame(res$prop) || is.matrix(res$prop))
})

# ----------------------
# CpG overlap warning/error tests
# ----------------------

test_that("run_houseman warns and proceeds if some CpGs are missing, errors if all are missing", {
  skip_if_not_installed("FlowSorted.Blood.EPIC")
  skip_if_not_installed("minfiData")
  skip_if_not_installed("minfi")
  sig <- methyldeconv::get_houseman_signature_matrix()
  subset_cpgs <- sig$CpGs[1:5]
  missing_cpgs <- paste0("cg_missing_", 1:5)
  methyl_set <- minfiData::MsetEx
  expect_warning(
    res <- methyldeconv::run_houseman(methyl_set, cpg_subset = c(subset_cpgs, missing_cpgs)),
    "not present in the Houseman signature matrix"
  )
  expect_error(
    methyldeconv::run_houseman(methyl_set, cpg_subset = missing_cpgs),
    "None of the specified cpg_subset CpGs are present"
  )
})

# ----------------------
# MethAtlas tests
# ----------------------

test_that("get_methatlas_signature_matrix returns a data frame with expected structure", {
  ref_path <- system.file("reference_atlas.csv", package = "methyldeconv")
  sig <- methyldeconv::get_methatlas_signature_matrix(ref_path)
  expect_s3_class(sig, "data.frame")
  expect_true(nrow(sig) > 0)
  expect_true(ncol(sig) > 0)
  expect_true("CpGs" %in% colnames(sig))
})

test_that("run_methatlas supports cpg_subset argument", {
  ref_path <- system.file("reference_atlas.csv", package = "methyldeconv")
  sig <- methyldeconv::get_methatlas_signature_matrix(ref_path)
  subset_cpgs <- sig$CpGs[1:10]
  beta <- matrix(runif(10 * 3), nrow = 10, ncol = 3)
  rownames(beta) <- subset_cpgs
  colnames(beta) <- paste0("Sample", 1:3)
  expect_error(
    methyldeconv::run_methatlas(beta, reference_atlas = ref_path, cpg_subset = subset_cpgs),
    NA
  )
})

# ----------------------
# CpG overlap warning/error tests
# ----------------------

test_that("run_methatlas warns and proceeds if some CpGs are missing, errors if all are missing", {
  ref_path <- system.file("reference_atlas.csv", package = "methyldeconv")
  sig <- methyldeconv::get_methatlas_signature_matrix(ref_path)
  subset_cpgs <- sig$CpGs[1:5]
  missing_cpgs <- paste0("cg_missing_", 1:5)
  beta <- matrix(runif(10 * 3), nrow = 10, ncol = 3)
  rownames(beta) <- c(subset_cpgs, missing_cpgs)
  colnames(beta) <- paste0("Sample", 1:3)
  expect_warning(
    methyldeconv::run_methatlas(beta, reference_atlas = ref_path, cpg_subset = c(subset_cpgs, missing_cpgs)),
    "not present in the reference atlas"
  )
  expect_error(
    methyldeconv::run_methatlas(beta, reference_atlas = ref_path, cpg_subset = missing_cpgs),
    "None of the specified cpg_subset CpGs are present"
  )
})

# ----------------------
# MethylResolver tests
# ----------------------

test_that("get_methylresolver_signature_matrix returns a data frame with CpGs column", {
  skip_if_not_installed("MethylResolver")
  sig <- methyldeconv::get_methylresolver_signature_matrix()
  expect_s3_class(sig, "data.frame")
  expect_true(nrow(sig) > 0)
  expect_true(ncol(sig) > 0)
  expect_true("CpGs" %in% colnames(sig))
})

test_that("run_methylresolver supports cpg_subset argument", {
  skip_if_not_installed("MethylResolver")
  sig <- methyldeconv::get_methylresolver_signature_matrix()
  n_celltypes <- ncol(sig) - 1 # exclude CpGs column
  n_cpgs <- max(25, 2 * n_celltypes)
  subset_cpgs <- sig$CpGs[1:n_cpgs]
  beta <- matrix(runif(n_cpgs * 3), nrow = n_cpgs, ncol = 3)
  rownames(beta) <- subset_cpgs
  colnames(beta) <- paste0("Sample", 1:3)
  res <- methyldeconv::run_methylresolver(beta, cpg_subset = subset_cpgs, alpha = 0.7)
  expect_true(is.list(res) || is.data.frame(res))
})

# ----------------------
# CpG overlap warning/error tests
# ----------------------

test_that("run_methylresolver warns and proceeds if some CpGs are missing, errors if all are missing", {
  skip_if_not_installed("MethylResolver")
  sig <- methyldeconv::get_methylresolver_signature_matrix()
  subset_cpgs <- sig$CpGs[1:5]
  missing_cpgs <- paste0("cg_missing_", 1:5)
  beta <- matrix(runif(10 * 3), nrow = 10, ncol = 3)
  rownames(beta) <- c(subset_cpgs, missing_cpgs)
  colnames(beta) <- paste0("Sample", 1:3)
  expect_warning(
    res <- methyldeconv::run_methylresolver(beta, cpg_subset = c(subset_cpgs, missing_cpgs), alpha = 0.7),
    "not present in the MethylResolver signature matrix"
  )
  expect_error(
    methyldeconv::run_methylresolver(beta, cpg_subset = missing_cpgs, alpha = 0.7),
    "None of the specified cpg_subset CpGs are present"
  )
})

# ----------------------
# MethylCC tests
# ----------------------

test_that("get_methylcc_signature_matrix returns a data frame with region columns", {
  skip_if_not_installed("FlowSorted.Blood.450k")
  skip_if_not_installed("methylCC")
  sig <- methyldeconv::get_methylcc_signature_matrix()
  expect_s3_class(sig, "data.frame")
  expect_true(nrow(sig) > 0)
  expect_true(ncol(sig) > 0)
  expect_true(all(c("seqnames", "start", "end") %in% colnames(sig)))
})

