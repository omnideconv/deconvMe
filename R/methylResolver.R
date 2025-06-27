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
#'
#' @export
#'
run_methylresolver <- function(beta_matrix, doPar = F, numCores = 1, alpha = seq(0.5,0.9,by = 0.05),
                               absolute = TRUE, purityModel = MethylResolver::RFmodel, seed = 1, cpg_subset = NULL){
  
  set.seed(seed)
  
  beta_matrix <- check_input_beta(beta_matrix)
  
  if (length(alpha) > 1){
    warning("MethylResolver may fail if multiple alpha values are provided. If this occurs, specify a single alpha value between 0.5 and 1.",
            immediate. = TRUE)
  }
  
  sig <- MethylResolver::MethylSig
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
  
  result_metrics <- result_methylresolver[,1:4]
  result_fractions <- result_methylresolver[,5:15]
  result_absolute <- result_methylresolver[,16:26]
  result_purity <- result_methylresolver[,27]
  
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