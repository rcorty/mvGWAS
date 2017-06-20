#' #' @title PermutationNullByDGSNP Function
#' #' @name PermutationNullByDGSNP
#' #' @author Robert Corty
#' #'
#' #' @description This function calculates the
#' #'
#' #' @param phenotype.name the phenotype to plot
#' #' @param snp.name the rsid of the SNP to stratify by
#' #' @param draw.ns should the plot show the number of individuals in each
#' #' genotype group?  Defaults to FALSE
#' #'
#' #' @return No return.  Just the plot.
#' #'
#' #'
#' PermutationNullByDGSNP <- function(phenotype.name,
#'                                    snp.name,
#'                                    file.names,
#'                                    box.cox.lambda = 1) {
#'
#'   snp <- GetDGSNPByName(file.names = file.names,
#'                         snp.name = snp.name)
#'
#'   phenos <- readRDS(file = file.names[['phen.file']])
#'
#'   perm.for.genos <- match(phenos$IID, names(snp))
#'
#'   snp <- snp[perm.for.genos]
#'
#'   df <- data_frame(snp = snp,
#'                    phen = if(box.cox.lambda == 0) {
#'                      log(phenos[[phenotype.name]])
#'                    } else {
#'                      phenos[[phenotype.name]]^box.cox.lambda
#'                    })
#'
#'
#'   recap <- df %>%
#'     group_by(snp) %>%
#'     summarise(count = n(),
#'               phen.mean = mean(phen, na.rm = TRUE),
#'               phen.sd = sd(phen, na.rm = TRUE),
#'               message = paste0('n = ', count, '\n',
#'                                'm = ', signif(phen.mean, 2), '\n',
#'                                's = ', signif(phen.sd, 2)))
#'
#'
#'
#' }