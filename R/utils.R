#' @title tryNA
#' @description silently get what you can out of an expression
#' @rdname internals
#'
#' @param expr the expression to try to evaluate
#'
#' @return If the expression evaluates without error, the evaluated expression.
#' If there is an error in evaluating the expression, NA.
#' Warnings are suppresssed.
#' @export
#'
tryNA <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NA,
                            finally = NA))
}


#' @title tryNULL
#' @description silently get what you can out of an expression
#' @rdname internals
#'
#'
#' @return If the expression evaluates without error, the evaluated expression.
#' If there is an error in evaluating the expression, NULL.
#' Warnings are suppresssed.
#' @export
#'
tryNULL <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NULL,
                            finally = NULL))
}





#' @title pull_GT
#' @description Pull the GT (maximum likelihood genotype) from a vcf file snp row.
#' @rdname internals
#'
#' @param snp_row the snp row from the VCF file (first cell is FORMAT)
#' @param min_gt_count minimum number of observations required to keep a genotype group in the study.
#'     Five is the default, but I advocate a slightly higher value for research (10 or 20 maybe).
#'
#' @return a data.frame where the first column is the ID and the second column is the genotype
#' @export
#'
pull_GT <- function(snp_row, min_gt_count) {

  if (names(snp_row)[1] != 'FORMAT') {
    stop('First column in genotype section of VCF file must be "FORMAT".')
  }

  snp_formats <- stringr::str_split(string = snp_row[1], pattern = ':')[[1]]
  gt_idx <- which(x = 'GT' == snp_formats)

  snp_row <- snp_row[-1]
  num_indiv <- length(snp_row)
  gts_raw <- stringr::str_split(string = snp_row, pattern = ':', simplify = TRUE)[,gt_idx]

  # if hets in both 'directions' are present, collapse them to one
  # todo: may need to deal with '0\1' or '0/1' at some point
  if (all('0|1' %in% gts_raw, '1|0' %in% gts_raw)) {
    gts <- replace(x = gts_raw, list = gts_raw == '1|0', values = '0|1')
  } else {
    gts <- gts_raw
  }

  # if there's any too-rare GT, drop it
  if (any(table(gts) < min_gt_count)) {
    bad_gts <- names(table(gts))[table(gts) < min_gt_count]
    gts <- replace(x = gts, list = which(gts == bad_gts), values = NA)
  }

  # if all the same GT, just return all 0's...
  # todo: check that this is kosher
  if (length(unique(gts)) == 1) {
    return(data.frame(ID = names(snp_row),
                      GT = rep(0, num_indiv),
                      stringsAsFactors = FALSE))
  } else {
    return(data.frame(ID = names(snp_row),
                      GT = factor(gts),
                      stringsAsFactors = FALSE))
  }

}



#' @title pull_DS
#' @description Pull the DS (alternate allele dosage) from a vcf file snp row.
#' @rdname internals
#'
#' @return a data.frame where the first column is the ID and the second column is the genotype
#' @export
#'
pull_DS <- function(snp_row) {

  if (names(snp_row)[1] != 'FORMAT') {
    stop('First column in genotype section of VCF file must be "FORMAT".')
  }

  snp_formats <- stringr::str_split(string = snp_row[1], pattern = ':')[[1]]
  ds_idx <- which(x = 'DS' == snp_formats)

  snp_row <- snp_row[-1]
  num_indiv <- length(snp_row)
  dosage_raw <- stringr::str_split(string = snp_row, pattern = ':', simplify = TRUE)[,ds_idx]
  dosage_num <- as.numeric(dosage_raw)

  return(data.frame(ID = names(snp_row),
                    DS = dosage_num,
                    stringsAsFactors = FALSE))
}



#' @title pull_GP
#' @description Pull the GP (genotype probabilities) from a vcf file snp row.
#' @rdname internals
#'
#' @return  a data.frame where the first column is the ID and the second and third columns are genotype probabilities
#' @export
#'
pull_GP <- function(snp_row) {

  if (names(snp_row)[1] != 'FORMAT') {
    stop('First column in genotype section of VCF file must be "FORMAT".')
  }

  snp_formats <- stringr::str_split(string = snp_row[1], pattern = ':')[[1]]
  gp_idx <- which(x = 'GP' == snp_formats)

  snp_row <- snp_row[-1]
  num_indiv <- length(snp_row)

  genoprob_string_vec <- stringr::str_split(string = snp_row, pattern = ':', simplify = TRUE)[,gp_idx]
  genoprobs_raw <- stringr::str_split(string = genoprob_string_vec, pattern = ',', simplify = TRUE)
  genoprobs_num <- apply(X = genoprobs_raw, MARGIN = 2, FUN = as.numeric)

  if (any(abs(apply(X = genoprobs_num, MARGIN = 1, FUN = sum) - 1) > 0.01)) {
    stop('Genoprobs dont add to 1.')
  }

  genoprob_df <- as.data.frame(genoprobs_num)
  names(genoprob_df) <- c('GP_ref', 'GP_het', 'GP_alt')
  genoprob_df$ID <- names(snp_row)

  return(genoprob_df[c('ID', 'GP_het', 'GP_alt')])
}
