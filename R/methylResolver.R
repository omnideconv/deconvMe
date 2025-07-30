#' Run MethylResolver
#'
#' @param beta_matrix a beta matrix with CpGs in rows and samples in columns
#' @param doPar Whether to use parallel processing to speed up the deconvolution computation if many samples are present. Default is FALSE. 
#' @param numCores Number of cores used for parallel processing to speed up the deconvolution of many samples. Requires doPar = TRUE. Default is 1. numCores = "auto" is max number of cores available minus one. 
#' @param alpha Set the alpha parameter for LTS deconvolution. This is the fraction of optimal CpGs from the signature matrix which are used for deconvolution. Must be between 0 and 1. Users can specify a vector or a single number. If a vector is specified, a grid search of the values is conducted and the alpha value that results in the lowest RMSE between the original and reconstructed mixture is selected. Default is seq(0.5,0.9,by = 0.05). 
#' @param absolute Whether to compute tumor purity and absolute cell type fractions. Default is TRUE. 
#' @param purityModel Random Forest model to predict mixture purity (unknown content) which allows the calculation of absolute cell type fractions. Required if absolute is TRUE. Default is our RF model trained on the consensus purity estimate (CPE) using TCGA data.  
#' @param seed fixed seed to account for RNG influences
#' @param cpg_subset Optional character vector of CpGs to subset the signature matrix. Default: NULL (use all CpGs in the signature).
#' @param reference Either the built-in MethylResolver signature matrix or a data.frame of reference CpGs with a 'CpGs' column and cell types as other columns.
#'
#' @export
#'
run_methylresolver <- function(beta_matrix, doPar = F, numCores = 1, alpha = seq(0.5,0.9,by = 0.05),
                               absolute = TRUE, purityModel = MethylResolver::RFmodel, seed = 1, cpg_subset = NULL, reference = NULL){
  
  set.seed(seed)
  
  beta_matrix <- check_input_beta(beta_matrix)
  
  if (length(alpha) > 1){
    warning("MethylResolver may fail if multiple alpha values are provided. If this occurs, specify a single alpha value between 0.5 and 1.",
            immediate. = TRUE)
  }
  
  # Handle reference parameter
  if (is.null(reference)) {
    sig <- MethylResolver::MethylSig
  } else if (is.data.frame(reference)) {
    if (!'CpGs' %in% colnames(reference)) {
      stop("No 'CpGs' column in custom reference data.frame.")
    }
    sig <- reference[, -which(colnames(reference) == 'CpGs')]
    rownames(sig) <- reference$CpGs
    colnames(sig) <- gsub(pattern = '-', replacement = '_', x = colnames(sig))
    
    warning('When a external signature matrix is provided, MethylResolver can no longer estimate tumor purity and absolute will be automatically set to FALSE.')
    absolute <- FALSE
  } else {
    stop("reference must be either NULL (use built-in) or a data.frame")
  }
  
  if (!is.null(cpg_subset)) {
    overlap <- intersect(cpg_subset, rownames(sig))
    missing <- setdiff(cpg_subset, rownames(sig))
    if (length(overlap) == 0) {
      stop("None of the specified cpg_subset CpGs are present in the MethylResolver signature matrix.")
    }
    if (length(missing) > 0) {
      msg <- paste0("Warning: The following CpGs are not present in the MethylResolver signature matrix and will be ignored: ",
        paste(head(missing, 10), collapse=", "))
      if (length(missing) > 10) msg <- paste0(msg, ", ...")
      warning(msg)
    }
    sig <- sig[rownames(sig) %in% overlap, , drop = FALSE]
  }
  
  result_methylresolver <- MethylResolver::MethylResolver(methylMix = beta_matrix, 
                                                          methylSig = sig, 
                                                          betaPrime = FALSE,
                                                          doPar = doPar, 
                                                          numCores = numCores, 
                                                          alpha = alpha, 
                                                          absolute = absolute, 
                                                          purityModel = purityModel)
  
  result_metrics <- result_methylresolver[,c('RMSE1', 'R1', 'RMSE2', 'R2')]
  result_fractions <- result_methylresolver[,colnames(sig)]
  if(absolute){
    result_absolute <- result_methylresolver[,paste0('abs_',colnames(sig))]
    result_purity <- result_methylresolver[['Purity']]
  }else{
    result_absolute <- NULL
    result_purity <- NULL
  }
  
  return(list(result_metrics=result_metrics,
              result_fractions=result_fractions,
              result_absolute=result_absolute,
              result_purity=result_purity))
}

#' Get MethylResolver Signature Matrix
#'
#' Returns the signature matrix used by MethylResolver.
#' @return Signature matrix as tibble with CpGs in rows (column 'CpGs') and cell types in columns
#' @export
get_methylresolver_signature_matrix <- function() {
  as.data.frame(MethylResolver::MethylSig) |> 
    tibble::rownames_to_column("CpGs")
}