library(deconvMe)

test_that("deconvolute errors on unsupported method", {
  expect_error(deconvolute(methyl_set = NULL, method = "notamethod"))
})

test_that("deconvolute errors on missing methyl_set", {
  expect_error(deconvolute(methyl_set = NULL, method = "epidish"))
})

test_that("deconvolute_combined errors on unsupported method", {
  # Use a dummy methyl_set, expect error for unsupported method
  expect_error(deconvolute_combined(methyl_set = NULL, methods = c("epidish", "notamethod")))
})

test_that("run_epidish errors on missing beta_matrix", {
  expect_error(run_epidish(beta_matrix = NULL))
})

test_that("run_houseman errors on missing methyl_set", {
  expect_error(run_houseman(methyl_set = NULL))
})

test_that("run_methylcc errors on missing methyl_set", {
  expect_error(run_methylcc(methyl_set = NULL))
})

test_that("run_methylresolver errors on missing beta_matrix", {
  expect_error(run_methylresolver(beta_matrix = NULL))
})

test_that("run_methatlas errors on missing beta_matrix", {
  expect_error(run_methatlas(beta_matrix = NULL))
}) 