library(methyldeconv)

test_that("results_boxplot runs without error on valid input", {
  df <- as.data.frame(matrix(runif(20), nrow=4))
  colnames(df) <- c("A","B","C","D","E")
  rownames(df) <- paste0("Sample", 1:4)
  expect_s3_class(results_boxplot(df), "ggplot")
})

test_that("results_barplot runs without error on valid input", {
  df <- as.data.frame(matrix(runif(20), nrow=4))
  colnames(df) <- c("A","B","C","D","E")
  rownames(df) <- paste0("Sample", 1:4)
  expect_s3_class(results_barplot(df), "ggplot")
})

test_that("results_aggregated_boxplot runs without error on valid input", {
  df <- data.frame(
    sample = rep(paste0("Sample", 1:4), 2),
    method = rep(c("epidish", "houseman"), each=4),
    celltype = rep(c("A","B"), each=2, times=2),
    value = runif(8)
  )
  expect_s3_class(results_aggregated_boxplot(df), "ggplot")
})

test_that("results_boxplot errors on invalid input", {
  expect_error(results_boxplot(data.frame(a=1:3, b=4:6)), NA)
  expect_error(results_boxplot(NULL))
})

test_that("results_barplot errors on invalid input", {
  expect_error(results_barplot(data.frame(a=1:3, b=4:6)), NA)
  expect_error(results_barplot(NULL))
})

test_that("results_aggregated_boxplot errors on invalid input", {
  expect_error(results_aggregated_boxplot(data.frame(a=1:3, b=4:6)))
  expect_error(results_aggregated_boxplot(NULL))
}) 