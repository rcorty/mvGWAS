#' @title CoValidatePhenotypesAndGenotypes
#' @name CoValidatePhenotypesAndGenotypes
#' @param phenTable, a table of phenotypes to be validated
#' @param genTable, a table of genotypes to be validated
#' @return TRUE if phenTable and genTable are co-valid for use in GWASviaDGLM
#' @author Robert Corty

CoValidatePhenotypesAndGenotypes <- function(phenTable, genTable) {


  return(list(phenTable = phenTable,
              genTable = genTable))
}
