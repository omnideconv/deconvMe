
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

#' Run selected set of methods and average results
#'
#' @param methyl_set A minfi MethylSet
#' @param array type of methylation array that was used. possible options are '450k' and 'EPIC'
#' @param methods list of methods (>1) that will be applied to the methyl set
#' @param scale_results   Whether the deconvolution results should be rescaled.
#'   Negative values will be set to 0, and the estimates will be normalized to sum to 1 per sample.
#'   Defaults to FALSE.
#' @param ... Additional parameters, passed to the algorithm used. See individual method documentations for details.
#'   
#' @return dataframe with results of all selected methods as well as the combined estimates
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
