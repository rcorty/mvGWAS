#' #' @title GenomeScan
#' #' @name GenomeScan
#' #' @param phenTable, a table of phenotypes to be validated
#' #' @param genTable, a table of genotypes to be validated
#' #' @return object of class c('tbl_df', 'data.frame') with []
#' #' @author Robert Corty
#'
#' GenomeScan <- function(meanFormula,
#'                        varFormula,
#'                        file.names,
#'                        verbose = FALSE) {
#'
#'   phenTable <- ValidatePhenotypes(phenTable = phenTable)
#'
#'   genTable <- ValidateGenotypes(genTable = genTable)
#'
#'   list[phenTable, genTable] <-
#'     CoValidatePhenotypesAndGenotypes(phenTable = validPhenTable,
#'                                      genTable = validGenTable)
#'
#'   num.snps <- ncol(genTable)
#'   result <- data_frame(mean.est = rep(NA, num.snps),
#'                        mean.se = rep(NA, num.snps),
#'                        var.est = rep(NA, num.snps),
#'                        var.se = rep(NA, num.snps),
#'                        n = rep(NA, num.snps))
#'   for (snp.idx in 1:num.snps) {
#'
#'     mappingTable <- cbind(phenTable,
#'                           genTable[[snp.idx]])
#'     names(mappingTable)[ncol(mappingTable)] <- 'snp'
#'
#'     fit <- dglm(formula = meanFormula,
#'                 dformula = varFormula,
#'                 data = mappingTable)
#'
#'     result[snp.idx, c('mean.est', 'mean.se')] <-
#'       summary(fit)$coef['snp', c('Estimate', 'Std. Error')]
#'     result[snp.idx, c('var.est', 'var.se')] <-
#'       summary(fit$dispersion.fit)$coef['snp', c('Estimate', 'Std. Error')]
#'     result[snp.idx, 'n'] <-
#'       sum(!is.na(phenTable[[meanFormula[[2]]]]) & !is.na(genTable[[snp.idx]]))
#'
#'   }
#'
#'   return(result)
#' }
