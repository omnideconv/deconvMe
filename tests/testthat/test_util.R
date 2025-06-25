library(methyldeconv)

test_that("normalize_deconv_results normalizes and sets negatives to zero", {
  mat <- matrix(c(0.2, -0.1, 0.9, 0.0), nrow=2, byrow=TRUE)
  colnames(mat) <- c("A", "B")
  res <- normalize_deconv_results(mat)
  expect_true(all(res >= 0))
  expect_equal(rowSums(res), c(1,1))
})

test_that("normalize_deconv_results handles single column", {
  mat <- matrix(c(0.2, 0.8), ncol=1)
  colnames(mat) <- "A"
  res <- normalize_deconv_results(mat)
  expect_equal(dim(res), c(2,1))
  expect_equal(colnames(res), "A")
})

test_that("rename_cell_types maps known and unknown names", {
  input <- c("CD8T", "CD4T", "Bcell", "NK", "Mono", "Neu", "foo")
  expected <- c("T cell CD8+", "T cell CD4+", "B cell", "NK cell", "Monocyte", "Neutrophil", "other")
  expect_equal(rename_cell_types(input), expected)
})

test_that("check_input_beta replaces NA with 0.5", {
  mat <- matrix(c(NA, 0.2, 0.3, NA), nrow=2)
  fixed <- suppressWarnings(check_input_beta(mat))
  expect_true(all(!is.na(fixed)))
  expect_true(all(fixed[is.na(mat)] == 0.5))
})

test_that("check_input errors on mismatched dimensions", {
  m1 <- matrix(1:4, nrow=2)
  m2 <- matrix(1:6, nrow=2)
  expect_error(check_input(m1, m2))
})

test_that("check_input errors on all NA", {
  m1 <- matrix(NA, nrow=2, ncol=2)
  expect_error(check_input(m1, m1))
})

test_that("check_input_mset errors on non-MethylSet", {
  expect_error(check_input_mset(matrix(1:4, nrow=2)))
}) 