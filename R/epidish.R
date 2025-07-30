#' run EpiDISH
#'
#' @param beta_matrix a beta matrix with CpGs in rows and samples in columns
#' @param mode Choice of a reference-based method ('RPC','CBS','CP')
#' @param reference Either a string (one of 'blood', 'breast', 'epithelial') or a data.frame of reference CpGs, i.e. representative sites, 
#' for a number of cell subtypes. If a matrix is provided, it will be used directly as the signature matrix.
#' @param maxit Only used in RPC mode, the limit of the number of IWLS iterations
#' @param nu.v Only used in CBS mode. It is a vector of several candidate nu values. nu is 
#' parameter needed for nu-classification, nu-regression, and 
#' one-classification in svm. The best estimation results among all candidate nu 
#' will be automatically returned.
#' @param constraint Only used in CP mode, you can choose either of 'inequality' or 'equality' 
#' normalization constraint. The default is 'inequality' (i.e sum of weights 
#' adds to a number less or equal than 1), which was implemented in 
#' Houseman et al (2012).
#' @param cpg_subset Optional character vector of CpGs to subset the signature matrix. Default: NULL (use all CpGs in the signature).
#'
#' @return CP-mode
#' A list with the following entries: estF: a matrix of the estimated fractions; 
#' ref: the reference centroid matrix used; dataREF: the subset of the input 
#' data matrix with only the probes defined in the reference matrix.
#' 
#' @return CBS-mode
#' A list with the following entries: estF: a matrix of the estimated fractions; 
#' nu: a vector of 'best' nu-parameter for each sample; 
#' ref: the reference centroid matrix used;
#' dataREF: the subset of the input data matrix with only the probes defined in the 
#' reference matrix.
#' 
#' @return RPC-mode
#' A list with the following entries: estF: a matrix of the estimated fractions;
#'  ref: the reference centroid matrix used; 
#' dataREF: the subset of the input data matrix with only the probes defined in the 
#' reference matrix.
#' @export
#'
run_epidish <- function(beta_matrix,
                        mode=c('RPC', 'CBS', 'CP'), 
                        reference=c('blood','breast','epithelial'), 
                        maxit = 50, nu.v = c(0.25, 0.5, 0.7), 
                        constraint = c("inequality", "equality"),
                        cpg_subset = NULL){
  
  beta_matrix <- check_input_beta(beta_matrix)

  if (length(mode) > 1) {
    mode <- mode[1]
    message(paste0(mode, " was chosen as default for \"mode\""))
  }
  # If reference is a character, use built-in, else use as matrix
  if (is.character(reference) && length(reference) > 1) {
    reference <- reference[1]
    message(paste0(reference, " was chosen as default for \"reference\""))
  }
  message(paste0("Starting EpiDISH deconvolution with mode ", mode, " ..."))

  if (is.character(reference)) {
    ref_mat <- switch(reference,
      'blood' = EpiDISH::centDHSbloodDMC.m,
      'breast' = EpiDISH::centEpiFibFatIC.m,
      'epithelial' = EpiDISH::centEpiFibIC.m,
      stop("Unknown reference type: ", reference)
    )
  } else if (is.data.frame(reference)) {
    
    # check if 'CpGs' column exists
    if(!'CpGs' %in% colnames(reference)){
      stop("No CpG column in custom reference.")
    }
    
    reference_cpgs <- reference[['CpGs']]
    ref_mat <- reference |> dplyr::select(-CpGs) |> as.matrix()
    rownames(ref_mat) <- reference_cpgs
    
  } else {
    stop("reference must be either a character string or a data.frame")
  }

  if (!is.null(cpg_subset)) {
    overlap <- intersect(cpg_subset, rownames(ref_mat))
    missing <- setdiff(cpg_subset, rownames(ref_mat))
    if (length(overlap) == 0) {
      stop("None of the specified cpg_subset CpGs are present in the reference matrix.")
    }
    if (length(missing) > 0) {
      msg <- paste0("Warning: The following CpGs are not present in the reference matrix and will be ignored: ",
        paste(head(missing, 10), collapse=", "))
      if (length(missing) > 10) msg <- paste0(msg, ", ...")
      warning(msg)
    }
    ref_mat <- ref_mat[rownames(ref_mat) %in% overlap, , drop = FALSE]
  }
  result_epidish <- EpiDISH::epidish(beta.m = beta_matrix,
                                     ref.m = ref_mat,
                                     method = mode,
                                     maxit = maxit,
                                     nu.v = nu.v,
                                     constraint = constraint)

  return(result_epidish)
}

#' Get EpiDISH Signature Matrix
#'
#' Returns the reference matrix used by EpiDISH for a given type.
#' @param reference One of 'blood', 'breast', 'epithelial'.
#' @return Signature matrix as tibble with CpGs in rows (column 'CpGs') and cell types in columns
#' @export
get_epidish_signature_matrix <- function(reference = c('blood', 'breast', 'epithelial')) {
  reference <- match.arg(reference)
  as.data.frame(switch(reference,
    'blood' = EpiDISH::centDHSbloodDMC.m,
    'breast' = EpiDISH::centEpiFibFatIC.m,
    'epithelial' = EpiDISH::centEpiFibIC.m
  )) |> tibble::rownames_to_column("CpGs")
}


