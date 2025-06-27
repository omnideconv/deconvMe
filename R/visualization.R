#' Function to plot deconvolution results as boxplots for each cell type
#'
#' @param result A data frame containing deconvolution results. Rows must correspond to samples, columns to cell types. Row names must be sample identifiers. The data frame can originate from `deconvolute()`, or from any external deconvolution method, as long as it matches this structure.
#'
#' @details
#' The input data frame must:
#'   - Have samples as rows and cell types as columns
#'   - Have row names set to sample identifiers
#'   - Contain only numeric values (cell type fractions or proportions)
#'
#' Example of a valid input:
#'   SampleID | CellTypeA | CellTypeB | CellTypeC
#'   -------- | --------- | --------- | ---------
#'   Sample1  |   0.2     |   0.5     |   0.3
#'   Sample2  |   0.1     |   0.7     |   0.2
#'
#' @export
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


#' Function to plot deconvolution results as barplots for each sample
#'
#' @param result A data frame containing deconvolution results. Rows must correspond to samples, columns to cell types. Row names must be sample identifiers. The data frame can originate from `deconvolute()`, or from any external deconvolution method, as long as it matches this structure.
#'
#' @details
#' The input data frame must:
#'   - Have samples as rows and cell types as columns
#'   - Have row names set to sample identifiers
#'   - Contain only numeric values (cell type fractions or proportions)
#'
#' Example of a valid input:
#'   SampleID | CellTypeA | CellTypeB | CellTypeC
#'   -------- | --------- | --------- | ---------
#'   Sample1  |   0.2     |   0.5     |   0.3
#'   Sample2  |   0.1     |   0.7     |   0.2
#'
#' @export
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

#' Function to plot aggregated deconvolution results as boxplots for each cell type and method
#'
#' @param result A data frame containing aggregated deconvolution results. Must have columns: 'sample', 'celltype', 'value', and 'method'. This can originate from `deconvolute_combined()` or any external source, as long as the structure matches.
#'
#' @details
#' The input data frame must:
#'   - Contain columns: 'sample', 'celltype', 'value', and 'method'
#'   - 'sample': sample identifier (character)
#'   - 'celltype': cell type name (character)
#'   - 'value': cell type fraction or proportion (numeric)
#'   - 'method': deconvolution method name (character)
#'
#' Example of a valid input:
#'   sample   | celltype | value | method
#'   -------- | -------- | ----- | --------
#'   Sample1  |   A      | 0.2   | method1
#'   Sample1  |   B      | 0.5   | method1
#'   Sample2  |   A      | 0.1   | method2
#'   Sample2  |   B      | 0.7   | method2
#'
#' @export
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



