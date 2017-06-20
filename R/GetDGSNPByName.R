#' #' @title GetDGSNPByName Function
#' #' @name GetDGSNPByName
#' #' @author Robert Corty
#' #'
#' #' @description This function accesses a snp by name based on a map.file,
#' #' a list of genotype.files, and a snp.name
#' #'
#' #' @param file.names a list of 3.
#' #'
#' #' @return the SNP in a named vector
#' #'
#' #' @importFrom magrittr "%>%"
#' #' @export
#' #'
#' GetDGSNPByName <- function(file.names, snp.name) {
#'
#'   map <- readRDS(file = file.names[['map.file']])
#'
#'   # figure out which geno.file needs to be loaded and load it
#'   chr <- map %>%
#'     dplyr::filter(rsID == snp.name) %>%
#'     dplyr::select(Chr) %>%
#'     as.character()
#'
#'   stopifnot(length(chr) == 1)
#'
#'   all.chrs <- sort(unique(map$Chr))
#'   chr.idx <- which(chr == all.chrs)
#'
#'   stopifnot(length(chr.idx) == 1)
#'
#'   genos <- readRDS(file = file.names[['geno.files']][chr.idx])
#'
#'   stopifnot('data.frame' %in% class(genos))
#'
#'   snp <- genos[, grep(pattern = snp.name, x = names(genos), fixed = TRUE)]
#'   names(snp) <- genos$IID
#'   attr(x = snp, which = 'chr') <- chr
#'
#'   return(snp)
#' }