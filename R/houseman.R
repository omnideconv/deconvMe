#' Run the improved Houseman method
#'
#' @param methyl_set A minfi MethylSet
#' @param array type of methylation array that was used. possible options are '450k' and 'EPIC'
#' @param compositeCellType Which composite cell type is being deconvoluted. Should be one of "Blood", "CordBloodCombined", "CordBlood", "CordBloodNorway", "CordTissueAndBlood", or "DLPFC". See details for preferred approaches.
#' @param processMethod Joint normalization/background correction for user and reference data. For MethylSet objects only "preprocessQuantile" is available. Set it to any minfi preprocessing function as a character if you want to override it, like "preprocessFunnorm"
#' @param probeSelect How should probes be selected to distinguish cell types? Options include: 1) "IDOL", (default) for using a customized set of probes obtained from IDOL optimization, available for Blood and Umbilical Cord Blood 2) "both", which selects an equal number (50) of probes (with F-stat p-value < 1E-8) with the greatest magnitude of effect from the hyper- and hypo-methylated sides, and 3) "any", which selects the 100 probes (with F-stat p-value < 1E-8) with the greatest magnitude of difference regardless of direction of effect. This according to minfi algorithm. Default input "auto" in minfi will use "any" for cord blood and "both" otherwise. Please see references for more details. 
#' @param cellTypes A vector of length K that contains the cell type names. Default: c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu"). Please notice that this library use Neutrophils instead of Granulocytes. See details for your library.
#' @param referencePlatform The platform for the reference dataset; options include c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC" (default), "IlluminaHumanMethylation27k"). If the input rgSet belongs to another platform, it will be converted using minfi function convertArray. Should not be changed by the user.
#' @param referenceset It is NULL by default. A custom reference RGChannelSet object (in quotes) if it is not a package installed. This option also allows the user to perform the deconvolution in closed computing clusters without internet access to ExperimentHub. For that download and save the reference and input the resulting object here. If using an installed reference package set to NULL. 
#' @param CustomCpGs a custom vector of probe names for cell deconvolution. For custom lists it should be a vector object (no quotes).
#' @param meanPlot Whether to plots the average DNA methylation across the cell-type discriminating probes within the mixed and sorted samples.
#' @param verbose Should the function be verbose?
#' @param lessThanOne Should the predictions be constrained to exactly one, in minfi default is FALSE, now you can select the option
#' @param cellCounts If cell counts are available (CBC, of flow sort) add a vector of lenght equal to the samples being deconvolved
#' @param cpg_subset Optional character vector of CpGs to subset the signature matrix. Default: NULL (use all CpGs in the signature).
#' @param ... Other arguments for preprocessQuantile or other normalizations
#'
#' @import FlowSorted.Blood.450k
#'
#' @export
#'
run_houseman <- function(methyl_set, array = c('450k','EPIC'),
                                compositeCellType=c('Blood','CordBloodCombined','CordBlood','CordBloodNorway','CordTissueAndBlood','DLPFC'),
                                processMethod = 'preprocessQuantile', probeSelect = c('IDOL','both','any'), cellTypes =c('CD8T','CD4T','NK','Bcell','Mono','Neu'),
                                referencePlatform = c('IlluminaHumanMethylationEPIC','IlluminaHumanMethylation450k','IlluminaHumanMethylation27k'),
                                referenceset = NULL, CustomCpGs = NULL, meanPlot = FALSE, verbose = TRUE, lessThanOne = FALSE, cellCounts = NULL, cpg_subset = NULL, ...){

  check_input_mset(methyl_set)
  options(matrixStats.useNames.NA = "deprecated")
  
  if (length(array) > 1) {
    array <- array[1]
    message(paste0(array, " was chosen as default for \"array\""))
  }
  if (length(compositeCellType) > 1) {
    compositeCellType <- compositeCellType[1]
    message(paste0(compositeCellType, " was chosen as default for \"compositeCellType\""))
  }
  if (length(referencePlatform) > 1) {
    referencePlatform <- referencePlatform[1]
    message(paste0(referencePlatform, " was chosen as default for \"referencePlatform\""))
  }
  if (length(probeSelect) > 1) {
    probeSelect <- probeSelect[1]
    message(paste0(probeSelect, " was chosen as default for \"probeSelect\""))
  }
  
  
  if(array == '450k'){
    minfi::`annotation`(methyl_set) <- c('array'='IlluminaHumanMethylation450k',
                                         'annotation'='ilmn12.hg19')
  }else if(array == 'EPIC'){
    minfi::`annotation`(methyl_set) <- c('array'='IlluminaHumanMethylationEPIC',
                                         'annotation'='ilm10b4.hg19')
  }

  # Use cpg_subset if provided, else CustomCpGs
  if (!is.null(cpg_subset)) {
    CustomCpGs <- cpg_subset
  }

  result_fsb <- FlowSorted.Blood.EPIC::estimateCellCounts2(rgSet = methyl_set,
                                                           compositeCellType = compositeCellType,
                                                           processMethod = processMethod, 
                                                           probeSelect = probeSelect, 
                                                           cellTypes = cellTypes, 
                                                           referencePlatform = referencePlatform, 
                                                           referenceset = referenceset, 
                                                           CustomCpGs = CustomCpGs, 
                                                           meanPlot = meanPlot, 
                                                           verbose = verbose, 
                                                           lessThanOne = lessThanOne, 
                                                           cellcounts = cellCounts, ...)
  
  return(result_fsb)
}

#' Get Houseman Signature Matrix
#'
#' Returns the IDOLOptimizedCpGs.compTable (signature matrix) used by Houseman for a given array type..
#' @return Signature matrix as tibble with CpGs in rows (column 'CpGs') and cell types in columns
#' @export
get_houseman_signature_matrix <- function() {
  as.data.frame(FlowSorted.Blood.EPIC::IDOLOptimizedCpGs.compTable) |> 
    tibble::rownames_to_column("CpGs")
}
