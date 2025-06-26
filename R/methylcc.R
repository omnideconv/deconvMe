#' run methylCC deconvolution
#'
#' @param methyl_set A minfi MethylSet
#' @param array type of methylation array that was used. possible options are '450k' and 'EPIC'
#' @param find_dmrs_object If the user would like to supply different differentially methylated regions, they can use the output from the find_dmrs function to supply different regions to estimatecc.
#' @param verbose TRUE/FALSE argument specifying if verbose messages should be returned or not. Default is TRUE.
#' @param epsilon Threshold for EM algorithm to check for convergence. Default is 0.01.
#' @param max_iter Maximum number of iterations for EM algorithm. Default is 100 iterations.
#' @param take_intersection TRUE/FALSE asking if only the CpGs included in object should be used to find DMRs. Default is FALSE.
#' @param include_cpgs TRUE/FALSE. Should individual CpGs be returned. Default is FALSE.
#' @param include_dmrs TRUE/FALSE. Should differentially methylated regions be returned. Default is TRUE.
#' @param init_param_method method to initialize parameter estimates. Choose between "random" (randomly sample) or "known_regions" (uses unmethyalted and methylated regions that were identified based on Reinus et al. (2012) cell sorted data.). Defaults to "random".
#' @param a0init Default NULL. Initial mean methylation level in unmethylated regions
#' @param a1init Default NULL. Initial mean methylation level in methylated regions
#' @param sig0init Default NULL. Initial var methylation level in unmethylated regions
#' @param sig1init Default NULL. Initial var methylation level in methylated regions
#' @param tauinit Default NULL. Initial var for measurement error
#' @param demo TRUE/FALSE. Should the function be used in demo mode to shorten examples in package. Defaults to FALSE.
#' @param seed fixed seed to account for RNG influences
#'
#' @import FlowSorted.Blood.450k  
#' @import FlowSorted.Blood.EPIC
#'
#' @return A object of the class estimatecc that contains information about the cell composition estimation (in the summary slot) and the cell composition estimates themselves (in the cell_counts slot).
#' @export
#'
run_methylcc <- function(methyl_set, array = c('450k','EPIC'),
                         find_dmrs_object = NULL, verbose = TRUE,
                         epsilon = 0.01, max_iter = 100, take_intersection = FALSE,
                         include_cpgs = FALSE, include_dmrs = TRUE,
                         init_param_method = "random", a0init = NULL, a1init = NULL,
                         sig0init = NULL, sig1init = NULL, tauinit = NULL, demo = FALSE,
                         seed = 1){
  set.seed(seed)
  options(matrixStats.useNames.NA = "deprecated")
  
  check_input_mset(methyl_set)
  
  if (length(array) > 1) {
    array <- array[1]
    message(paste0(array, " was chosen as default for \"array\""))
  }
  
  if(array == '450k'){
    minfi::`annotation`(methyl_set) <- c('array'='IlluminaHumanMethylation450k',
                                         'annotation'='ilmn12.hg19')
  }else if(array == 'EPIC'){
    minfi::`annotation`(methyl_set) <- c('array'='IlluminaHumanMethylationEPIC',
                                         'annotation'='ilm10b4.hg19')
  }
  genomic_methyl_set <- minfi::mapToGenome(methyl_set)
  
  cell_composition <-
    methylCC::estimatecc(
      object = genomic_methyl_set,
      find_dmrs_object = find_dmrs_object,
      verbose = verbose,
      epsilon = epsilon,
      max_iter = max_iter,
      take_intersection = take_intersection,
      include_cpgs = include_cpgs,
      include_dmrs = include_dmrs,
      init_param_method = init_param_method,
      a0init = a0init,
      a1init = a1init,
      sig0init = sig0init,
      sig1init = sig1init,
      tauinit = tauinit,
      demo = demo
    )
  result_methylCC <- methylCC::cell_counts(object = cell_composition)
  
  return(result_methylCC)
}

#' Get MethylCC Signature Matrix (DMRs)
#'
#' Runs find_dmrs on the default reference and returns the DMRs.
#' @return DMRs as dataframe (as used by methylCC)
#' @export
get_methylcc_signature_matrix <- function() {
  dmr_list <- methylCC:::.find_dmrs()
  signature_dmrs <- dmr_list$regions_all

  data.frame(signature_dmrs)

}
