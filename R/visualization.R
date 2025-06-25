#' Function to plot the results from `deconvolute()`  as boxplots for each celltype
#'
#' @param result result from `deconvolute()` 
#'
#' @export
#'
#' 
results_boxplot <- function(result){
  if (!is.data.frame(result)) {
    stop("Input to results_boxplot must be a data frame.")
  }
  if (ncol(result) < 1) {
    stop("Input data frame must have at least one column.")
  }
  if (is.null(rownames(result))) {
    stop("Input data frame must have row names.")
  }
  result |> 
    tibble::rownames_to_column(var = 'sample') |>
    tidyr::pivot_longer(cols = -sample, names_to = 'celltype') |>
    ggplot2::ggplot(mapping = ggplot2::aes(x=value, y= celltype, fill=celltype))+
      ggplot2::geom_boxplot()+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = 'top')
}


#' Function to plot the results from `deconvolute()`  as barplots for each sample
#'
#' @param result result from `deconvolute()` 
#'
#' @export
#'
#' 
results_barplot <- function(result){
  if (!is.data.frame(result)) {
    stop("Input to results_barplot must be a data frame.")
  }
  if (ncol(result) < 1) {
    stop("Input data frame must have at least one column.")
  }
  if (is.null(rownames(result))) {
    stop("Input data frame must have row names.")
  }
  result |> 
    tibble::rownames_to_column(var = 'sample') |>
    tidyr::pivot_longer(cols = -sample, names_to = 'celltype') |>
    ggplot2::ggplot(mapping = ggplot2::aes(x=value, y=sample, fill=celltype))+
      ggplot2::geom_col()+
      ggplot2::theme_bw()+
      ggplot2::theme(legend.position = 'top')
}

#' Function to plot the results from `deconvolute_combined()` as boxplots for each celltype and celltype
#'
#' @param result result from `deconvolute_combined()` 
#'
#' @export
#'
#' 
results_aggregated_boxplot <- function(result){
  if (!is.data.frame(result)) {
    stop("Input to results_aggregated_boxplot must be a data frame.")
  }
  required_cols <- c("celltype", "value", "method")
  missing_cols <- setdiff(required_cols, colnames(result))
  if (length(missing_cols) > 0) {
    stop(paste("Input data frame is missing required columns:", paste(missing_cols, collapse=", ")))
  }
  ggplot2::ggplot(result, mapping = ggplot2::aes(x=value, y= celltype, fill=method))+
    ggplot2::geom_boxplot()+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = 'top')
}



