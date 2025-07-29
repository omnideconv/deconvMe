
#' Check input beta matrix
#'
#' @param beta_matrix beta matrix
#'
#' @return same matrix, but fixed
#' @export
#'
check_input_beta <- function(beta_matrix){
  beta_matrix <- as.matrix(beta_matrix)
  
  if(any(is.na(beta_matrix))){
    n_na <- length(which(is.na(beta_matrix)))
    warning(paste0(n_na, ' NA values detected in your beta matrix. Replacing them with 0.5.'))
    beta_matrix[which(is.na(beta_matrix))] <- 0.5
  }
  
  return(beta_matrix)
}

#' Check input matrices
#'
#' @param meth matrix with quantities of methylated sites
#' @param unmeth matrix with quantities of methylated sites
#'
check_input <- function(meth, unmeth){
  if(!all(dim(meth) == dim(unmeth))){
    stop('Meth and Unmeth matrices have different dimensions, stopping.')
  }
  if(all(is.na(meth))){
    stop('All entries in Meth matrix are NA, stopping.')
  }
  if(all(is.na(unmeth))){
    stop('All entries in Unmeth matrix are NA, stopping.')
  }
}

#' Check input methylSet
#'
#' @param methyl_set minfi methylSet
#'
#'
check_input_mset <- function(methyl_set){
  if(!is(methyl_set, "MethylSet")){
    stop('Input is not a MethylSet, stopping.')
  }
  if(all(is.na(methyl_set))){
    stop('All entries in the MethylSet are NA, stopping.')
  }
}


#' Create beta matrix from a methylset matrices
#'
#' @param methyl_set minfi methylSet
#'
#' @return beta matrix
#'
create_beta <- function(methyl_set){
  ratio_set <- minfi::ratioConvert(methyl_set)
  
  return(minfi::getBeta(ratio_set))
}

create_genomicMethylSet <- function(beta_matrix, array_type){
  
  genomic_ratio_set <- makeGenomicRatioSetFromMatrix(
    mat = beta_matrix,
    pData = NULL
  )
  
  if(array_type == '450k'){
    minfi::`annotation`(genomic_ratio_set) <- c('array'='IlluminaHumanMethylation450k',
                                                'annotation'='ilmn12.hg19')
  }else if(array_type == 'EPIC'){
    minfi::`annotation`(genomic_ratio_set) <- c('array'='IlluminaHumanMethylationEPIC',
                                                'annotation'='ilm10b4.hg19')
  }
  
  genomic_methyl_set <- as(genomic_ratio_set, "GenomicMethylSet")
}


#' Normalize deconvolution result
#'
#' @param deconv_result The original deconvolution result
#'
#' @return A matrix with the rowsums of one and no negative values
#' @export
normalize_deconv_results <- function(deconv_result) {
  celltypes <- colnames(deconv_result)
  deconv_result[deconv_result < 0] <- 0
  deconv_result <- t(apply(deconv_result, 1, function(row) row / sum(row)))
  # Apply returns a vector when only supplied with one celltype. To counter it and return a matrix
  # and not a vector, this operation is needed
  if (length(celltypes) == 1) {
    deconv_result <- t(deconv_result)
    colnames(deconv_result) <- celltypes
  }
  return(deconv_result)
}


#' Rename Cell Types in Deconvolution Results
#'
#' Renames a set of cell type names to standardized names.
#' It uses predefined mappings to replace specific cell type names with more 
#' consistent labels (e.g., "CD8T" to "T cell CD8+"). Unrecognized names are labeled as "other".
#'
#' @param input_celltypes list of cell types
#' @return a new list of cell types with controlled vocabulary
#' @export
#' 
rename_cell_types <- function(input_celltypes){
  output_celltypes <- dplyr::recode(input_celltypes,
                                     "CD8T" = "T cell CD8+", 
                                     "CD8" = "T cell CD8+", 
                                     "CD8T-cells_EPIC" = "T cell CD8+", 
                                     
                                     "CD4T" = "T cell CD4+", 
                                     "CD4T-cells_EPIC" = "T cell CD4+", 
                                     "T cell CD4+" = "T cell CD4+",
                                     
                                     "B" = "B cell", 
                                     "Bcell" = "B cell", 
                                     "B-cells_EPIC" = "B cell",
                                     
                                     "NK" = "NK cell",
                                     "NK-cells_EPIC" = "NK cell",
                                     
                                     "Mono" = "Monocyte", 
                                     "Mon" = "Monocyte",
                                     "Monocytes_EPIC" = "Monocyte",
                                     
                                     "Neu" = "Neutrophil",
                                     "Neutro" = "Neutrophil",
                                     "Neutrophils" = "Neutrophil",
                                     "Neutrophils_EPIC" = "Neutrophil",
                                           
                                     .default = "other"
  )
  
  return(output_celltypes)
}


#' Initialize Python Environment with Miniconda
#'
#' Sets up a Miniconda environment with Python 3.8 for the 'r-deconvMe' package.
#' Installs Miniconda if not present, creates the 'r-deconvMe' Conda environment 
#' if it doesn't exist, and installs required Python packages (`numpy`, `pandas`, 
#' `scipy`, `matplotlib`).
#'
#' @details 
#' - Installs Miniconda if missing.
#' - Creates 'r-deconvMe' environment with Python 3.8 and installs dependencies.
#' - Activates the environment and displays Python configuration.
#'
#' @import reticulate
#' @export
#' 
#' @examples
#' \dontrun{
#'   init_python()
#' }
init_python <- function(){
  
  if (!dir.exists(reticulate::miniconda_path())) {
    message("Setting python version in miniconda to be 3.8")
    Sys.setenv(RETICULATE_MINICONDA_PYTHON_VERSION = 3.8)
    message("Setting up miniconda environment..")
    suppressMessages(reticulate::install_miniconda())
  }
  
  if (!(reticulate::condaenv_exists("r-deconvMe"))) {
    message("Create conda evironment 'r-deconvMe' for MethAtlas...")
    reticulate::conda_create("r-deconvMe", python_version = "3.8")
    message("Install all python dependencies...")
    reticulate::py_install(packages = c("numpy", "pandas", "scipy", "matplotlib") , envname = "r-deconvMe",  method = "conda", pip = T)
  }

  
  reticulate::use_condaenv(condaenv = "r-deconvMe", required = FALSE)
  reticulate::py_config()
  
}