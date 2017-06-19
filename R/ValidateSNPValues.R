#' @title ValidateSNPValues
#' @name ValidateSNPValues
#' @author Robert W. Corty
#'
#' @param snp.values the snp values
#' @param phen.df the phenotype data.frame
#'
#' @return nothing if the snp values are valid, errors if they're not
#' @export
#'
ValidateSNPValues <- function(snp.values, phen.df) {

  if (any(is.na(match(names(snp.values), phen.df$IID)))) {
    stop('mismtach between names of snp.values and IID column in phenotype df')
  }
  if (any(is.na(match(phen.df$IID, names(snp.values))))) {
    stop('mismtach between names of snp.values and IID column in phenotype df')
  }
}