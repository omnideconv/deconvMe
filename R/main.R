
#' List of supported deconvolution methods
#'
#' The methods currently supported are
#' `EpiDISH`, `Houseman`, `MethylCC`, `MethylResolver`, `MethAtlas`
#'
#' The object is a named vector. The names correspond to the display name of the method,
#' the values to the internal name.
#'
#' @export
deconvolution_methods <- c(
  "EpiDISH" = "epidish", "Houseman" = "houseman", "MethylCC" = "methylcc", 
  "MethylResolver" = "methylresolver", "MethAtlas" = "methatlas"
)


#' Deconvolution with deconvMe
#'
#' @param methyl_set A minfi MethylSet
#' @param method A string specifying the method. Supported methods are 'epidish', 'houseman', 'methylcc', 'methylresolver', 'methatlas'
#' @param scale_results   Whether the deconvolution results should be rescaled.
#'   Negative values will be set to 0, and the estimates will be normalized to sum to 1 per sample.
#'   Defaults to FALSE.
#' @param ... Additional parameters, passed to the algorithm used. See individual method documentations for details.
#'
#' @return A matrix with the probabilities of each cell-type for each individual. Rows are
#' individuals, columns are cell types.
#' @export
#'
#' @examples 
#' 
#' ex_data <- minfiData::MsetEx
#' 
#' result <- deconvolute(ex_data, method='epidish')
deconvolute <- function(methyl_set, method=deconvolution_methods, scale_results = FALSE, ...){
  
  options(matrixStats.useNames.NA = "deprecated")
  
  if (length(method) > 1) {
    stop(
      "Please only specify one method and not ", length(method), ": ",
      paste(method, collapse = ", ")
    )
  }
  
  if (method %in% names(deconvolution_methods)) {
    method <- deconvolution_methods[[method]]
  }
  
  check_input_mset(methyl_set)
  beta_matrix <- create_beta(methyl_set)
  
  method <- tolower(method)
  
  result <- switch (method,
    epidish = run_epidish(beta_matrix, ...)$estF,
    houseman = run_houseman(methyl_set, ...)$prop,
    methylcc = as.matrix(run_methylcc(methyl_set, ...)),
    methylresolver = as.matrix(run_methylresolver(beta_matrix, ...)$result_fractions),
    methatlas = run_methatlas(beta_matrix, ...)
  )
  
  if(!is.null(result)){
    # Scale the results to sum up to 1
    if (scale_results) {
      deconv <- normalize_deconv_results(result)
    }
    # Alphabetical order of celltypes
    result <- result[, order(colnames(result)), drop = FALSE]
    
    # transform to dataframe
    result <- result |> data.frame(check.names = F)
  }
  
  return(result)
}

#' Run multiple deconvolution methods and create aggregated results
#'
#' This function runs multiple cell-type deconvolution methods on the same methylation data and creates both individual method results and an aggregated (averaged) estimate. The aggregation approach can help reduce method-specific biases and provide more robust cell-type proportion estimates.
#'
#' \strong{How it works:}
#' \itemize{
#' \item Runs each specified method independently on the methylation data
#' \item Standardizes cell-type names across methods using \code{rename_cell_types()}
#' \item For MethylResolver specifically, combines Tnaive and Tmem into "T cell CD4+" to match other methods
#' \item Calculates the mean estimate for each cell type across all methods (aggregated results)
#' \item Returns both individual method results and aggregated results in a long-format data frame
#' }
#'
#' \strong{Cell-type standardization:}
#' The function uses \code{rename_cell_types()} to standardize cell-type names across different methods. This mapping includes:
#' \itemize{
#' \item CD8T/CD8/CD8T-cells_EPIC → "T cell CD8+"
#' \item CD4T/CD4T-cells_EPIC → "T cell CD4+"
#' \item B/Bcell/B-cells_EPIC → "B cell"
#' \item NK/NK-cells_EPIC → "NK cell"
#' \item Mono/Mon/Monocytes_EPIC → "Monocyte"
#' \item Neu/Neutro/Neutrophils/Neutrophils_EPIC → "Neutrophil"
#' \item Any unrecognized cell types → "other"
#' }
#'
#' \strong{Meaning of 'other' cell types:}
#' The "other" category includes:
#' \itemize{
#' \item Cell types that are not in the standardized mapping above
#' \item Method-specific cell types that don't have direct equivalents in other methods
#' \item Rare or specialized cell populations that are only detected by certain methods
#' }
#'
#' \strong{Limitations of the aggregation approach:}
#' \itemize{
#' \item \strong{Method heterogeneity:} Different methods use different algorithms, reference datasets, and cell-type definitions, which may not be directly comparable
#' \item \strong{Missing data handling:} If a method fails to estimate a particular cell type, it may be excluded from the aggregation, potentially biasing results
#' \item \strong{Equal weighting:} The current implementation gives equal weight to all methods, regardless of their individual performance or reliability
#' \item \strong{Cell-type coverage:} Methods may detect different sets of cell types, leading to incomplete coverage in aggregated results
#' \item \strong{No confidence intervals:} The aggregation provides point estimates without uncertainty quantification
#' }
#'
#' \strong{Best practices:}
#' \itemize{
#' \item Use methods that are well-validated for your specific tissue type and experimental design
#' \item Consider the biological context when interpreting aggregated results
#' \item Validate results against independent measurements when possible
#' \item Be cautious when aggregating methods with very different cell-type definitions
#' \item Consider using the individual method results to assess consistency across methods
#' }
#'
#' @param methyl_set A minfi MethylSet
#' @param array type of methylation array that was used. possible options are '450k' and 'EPIC'
#' @param methods list of methods (>1) that will be applied to the methyl set
#' @param scale_results   Whether the deconvolution results should be rescaled.
#'   Negative values will be set to 0, and the estimates will be normalized to sum to 1 per sample.
#'   Defaults to FALSE.
#' @param ... Additional parameters, passed to the algorithm used. See individual method documentations for details.
#'   
#' @return A data frame with columns: sample, method, celltype, value. Contains results from all individual methods plus an 'aggregated' method that averages the estimates across methods.
#' @export
#' @examples 
#' 
#' ex_data <- minfiData::MsetEx
#' 
#' result <- deconvolute_combined(ex_data, methods=c('epidish','houseman'))
deconvolute_combined <- function(methyl_set, array = c('450k','EPIC'), methods, scale_results = FALSE, ...){
  
  if(any(!methods %in% deconvolution_methods)){
    stop(paste0('At least one of your selected methods is not supported by deconvMe. Please check your spelling, supported methods are: ',
                'epidish, houseman, methylcc, methylresolver, methatlas'))
  }
  
  result <- lapply(methods, function(m){
    df <- deconvolute(methyl_set = methyl_set, method = m, scale_results = scale_results, ...) |>
      tibble::rownames_to_column(var = 'sample')
    
    # for methylresolver, combine Tnaive and Tmem to CD4+ cells
    if(m == 'methylresolver'){
      df <- cbind(df, "T cell CD4+" = df[, "Tmem"] + df[, "Tnaive"])
    }

    return(df)
  })
  
  names(result) <- methods
  method_results <- dplyr::bind_rows(result, .id = 'method') |>
    tidyr::pivot_longer(cols = -c(sample, method), names_to = 'celltype') |>
    dplyr::mutate(celltype = rename_cell_types(celltype)) |>
    dplyr::filter(!is.na(value)) 
  
  
  combined_results <- method_results |>
    dplyr::group_by(sample, celltype) |>
    dplyr::summarize(value = mean(value), .groups = 'drop') |>
    dplyr::mutate(method = 'aggregated')
  
  full_results <- dplyr::bind_rows(method_results, combined_results)

  return(full_results)
}
