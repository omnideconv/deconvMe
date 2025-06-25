#' run MethAtlas
#'
#' @param beta_matrix a beta matrix with CpGs in rows and samples in columns
#' @param reference_atlas Path to a csv file that saves a reference matrix with CpGs as rows and cell types as columns.
#'                        The default (tissue-wide) reference file is stored in 'inst/reference_atlas.csv'. 
#' @param temp_dir Path to directory where the beta matrix will be saved as a csv file.
#' @param out_dir Path to output directory. Output will be a csv file and a png representing the cell type fractions.
#' @param use_epic_reference The MethAtlas has a whole-tissue reference or a immunecell-specific reference that is optimized for EPIC arrays (which is a subset of the whole-tissue reference)
#' 
#' @export
#'
run_methatlas <- function(beta_matrix, reference_atlas = system.file("reference_atlas.csv", package = "methyldeconv"), temp_dir = NULL, out_dir = NULL, use_epic_reference=FALSE){
  # check if python is installed, else install
  init_python()
  
  # set up temporary nd output directories
  tmp_dir <- temp_dir
  if (is.null(temp_dir)) {
    tmp_dir <- tempdir()
    dir.create(tmp_dir, showWarnings = FALSE)
  }
  
  if (is.null(out_dir)) {
    out_dir <- tmp_dir
  } 
  
  if (!dir.exists(out_dir)){
    dir.create(out_dir, showWarnings = FALSE)
  }
  
  # create a beta matrix from the methyl_set and save to temporary folder
  beta_matrix <- check_input_beta(beta_matrix)
  beta_path = paste0(tmp_dir, "/beta.csv")
  write.csv(beta_matrix, beta_path)
  
  # subset reference if applicable
  if(use_epic_reference){
    reference_atlas <- system.file("reference_atlas_epic.csv", package = "methyldeconv")
  }
  
  # run meth_atlas
  system(paste("python", system.file("deconvolve.py", package = "methyldeconv")," -a", reference_atlas, beta_path, "--out", out_dir))
  
  # read the results to provide as data frame
  as.matrix(t(read.csv(paste0(out_dir, "/beta_deconv_output.csv"),
                      row.names = 1, check.names = FALSE
  )))
}

#' Get MethAtlas Signature Matrix
#'
#' Returns the reference matrix used by MethAtlas from a CSV file.
#' @param reference_atlas Path to a csv file (default: package reference_atlas.csv)
#' @return Signature matrix as tibble with CpGs in rows (column 'CpGs') and cell types in columns
#' @export
get_methatlas_signature_matrix <- function(reference_atlas = system.file("reference_atlas.csv", package = "methyldeconv")) {
  read.csv(reference_atlas, check.names = FALSE)
}
  