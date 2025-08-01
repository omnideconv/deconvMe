#' run MethAtlas
#'
#' @param beta_matrix a beta matrix with CpGs in rows and samples in columns
#' @param reference Either a path to a csv file that saves a reference matrix with CpGs as rows and cell types as columns,
#'                  or a data.frame of reference CpGs with a 'CpGs' column and cell types as other columns.
#'                  The default (tissue-wide) reference file is stored in 'inst/reference_atlas.csv'. 
#' @param temp_dir Path to directory where the beta matrix will be saved as a csv file.
#' @param out_dir Path to output directory. Output will be a csv file and a png representing the cell type fractions.
#' @param use_epic_reference The MethAtlas has a whole-tissue reference or a immunecell-specific reference that is optimized for EPIC arrays (which is a subset of the whole-tissue reference)
#' @param cpg_subset Optional character vector of CpGs to subset the signature matrix. Default: NULL (use all CpGs in the signature).
#' 
#' @export
#'
run_methatlas <- function(beta_matrix, reference = system.file("reference_atlas.csv", package = "deconvMe"), temp_dir = NULL, out_dir = NULL, use_epic_reference=FALSE, cpg_subset = NULL){
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
  
  # Handle reference parameter
  if (is.character(reference)) {
    # If it's a character, treat as file path and check if it exists
    if(!file.exists(reference)) {
      stop("Reference file does not exist: ", reference)
    }
    
    # If it's epic reference, use the epic reference file
    if(use_epic_reference){
      reference <- system.file("reference_atlas_epic.csv", package = "deconvMe")
    }
    ref_df <- read.csv(reference, check.names = FALSE)
  } else if (is.data.frame(reference)) {
    # If it's a data.frame, use it directly
    if (!'CpGs' %in% colnames(reference)) {
      stop("No 'CpGs' column in custom reference data.frame.")
    }
    ref_df <- reference
  } else {
    stop("reference must be either a character string (file path) or a data.frame")
  }
  
  if (!is.null(cpg_subset)) {
    overlap <- intersect(cpg_subset, ref_df$CpGs)
    missing <- setdiff(cpg_subset, ref_df$CpGs)
    if (length(overlap) == 0) {
      stop("None of the specified cpg_subset CpGs are present in the reference atlas.")
    }
    if (length(missing) > 0) {
      msg <- paste0("Warning: The following CpGs are not present in the reference atlas and will be ignored: ",
        paste(head(missing, 10), collapse=", "))
      if (length(missing) > 10) msg <- paste0(msg, ", ...")
      warning(msg)
    }
    ref_df <- ref_df[ref_df$CpGs %in% overlap, , drop = FALSE]
    # Write the subsetted reference to a temp file
    reference_path <- paste0(tmp_dir, "/reference_atlas_subset.csv")
    write.csv(ref_df, reference_path, row.names = FALSE)
  }else{
    # Write the reference to a temp file
    reference_path <- paste0(tmp_dir, "/reference_atlas_temp.csv")
    write.csv(ref_df, reference_path, row.names = FALSE)
  }
  
  # run meth_atlas
  system(paste("python", system.file("deconvolve.py", package = "deconvMe")," -a", reference_path, beta_path, "--out", out_dir))
  
  # read the results to provide as data frame
  as.matrix(t(read.csv(paste0(out_dir, "/beta_deconv_output.csv"),
                      row.names = 1, check.names = FALSE
  )))
}

#' Get MethAtlas Signature Matrix
#'
#' Returns the reference matrix used by MethAtlas from a CSV file.
#'
#' @return Signature matrix as tibble with CpGs in rows (column 'CpGs') and cell types in columns
#' @export
get_methatlas_signature_matrix <- function() {
  reference_atlas = system.file("reference_atlas.csv", package = "deconvMe")
  unique(read.csv(reference_atlas, check.names = FALSE))
}
  