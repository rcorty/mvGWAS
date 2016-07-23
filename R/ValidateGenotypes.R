#' @title ValidateGenotypes
#' @name ValidateGenotypes
#' @param genTable, a table of genotypes to be validated
#' @return TRUE if genTable is a valid genotype table for GWASviaDGLM, FALSE otherwise
#' @author Robert Corty

ValidateGenotypes <- function(genTable) {

  genoVals <- genTable %>% select(-id)

  if (min(genoVals) < 0) { stop('Genotype less than 0 observed. Genotypes must be in [0, 2].')}
  if (max(genoVals) > 2) { stop('Genotype greater than 2 observed. Genotypes must be in [0, 2].')}

  return(genTable)
}
