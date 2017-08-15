#' @title DGLM-norm
#'
#' @param m.form mean formula
#' @param d.form variance formula
#' @param data data
#' @param maxiter max number of iterations to do
#' @param conv converges when LL changes by less than x
#'
#' @description this is a function that needs to be documented
#'
#' @return a list of length 3
#' @export
#'
DGLM_norm <- function(m.form, d.form, data, maxiter=20, conv=1e-6) {
  X.mean <- model.matrix(m.form, data = data)
  X.disp <- model.matrix(d.form, data = data)
  y.name <- all.vars(m.form)[1]
  y <- data[[y.name]]
  w <- rep(1, nrow(data))
  convergence <- 1
  iter <- 0
  while (convergence > conv & iter < maxiter) {
    iter <- iter + 1
    w.old <- w
    glm1 <- lm(y~.-1, weights = w, data=data.frame(X.mean))
    res <- resid(glm1)
    q <- hatvalues(glm1)
    y2 <- res^2
    glm2 <- glm(y2~.-1, family = stats::Gamma(link = log), weights = (1 - q)/2, data = data.frame(X.disp))
    w <- 1/fitted(glm2)
    convergence <- (max(abs(w.old - w)) + (summary(glm1)$sigma - 1) )
  }
  return(list(mean = glm1, disp = glm2, iter = iter))
}


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
    gts <- replace(x = gts, list = which(gts %in% bad_gts), values = NA)
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





#' fit_dglm
#'
#' @param mean_alt_formula yada
#' @param var_alt_formula yada
#' @param mean_null_formula yada
#' @param var_null_formula yada
#' @param this_locus_df yada
#'
#' @return a one-row data_frame
#' @export
#'
fit_dglm <- function(mean_alt_formula,
                     var_alt_formula,
                     mean_null_formula,
                     var_null_formula,
                     this_locus_df) {

  alt_fit <- tryNULL(DGLM_norm(m.form = mean_alt_formula,
                               d.form = var_alt_formula,
                               data = this_locus_df))

  # if we couldn't fit the alt model, move on
  if (is.null(alt_fit)) {
    return(dplyr::data_frame(LR_mean = NA,
                             df_mean = NA,
                             LR_var = NA,
                             df_var = NA,
                             LR_joint = NA,
                             df_joint = NA))
  }

  if ('DS' %in% all.vars(mean_alt_formula)) {

    beta_DS_mean <- tryNA(coef(summary(alt_fit$mean))['DS', 'Estimate'])
    se_DS_mean <- tryNA(coef(summary(alt_fit$mean))['DS', 'Std. Error'])
  }

  if ('DS' %in% all.vars(var_alt_formula)) {

    beta_DS_var <- tryNA(coef(summary(alt_fit$dispersion.fit))['DS', 'Estimate'])
    se_DS_var <- tryNA(coef(summary(alt_fit$dispersion.fit))['DS', 'Std. Error'])
  }


  # mean test
  if (exists(x = 'mean_null_formula')) {

    mean_null_fit <- tryNULL(DGLM_norm(m.form = mean_null_formula,
                                       d.form = var_alt_formula,
                                       data = this_locus_df))

    if (is.null(mean_null_fit)) {
      LR_mean <- df_mean <- NA
    } else {
      LR_mean <- as.numeric(2*logLik(alt_fit$mean) - 2*logLik(mean_null_fit$mean))
      df_mean <- alt_fit$mean$rank - mean_null_fit$mean$rank
    }
  }


  # var test
  if (exists(x = 'var_null_formula')) {

    var_null_fit <- tryNULL(DGLM_norm(m.form = mean_alt_formula,
                                      d.form = var_null_formula,
                                      data = this_locus_df))

    if (is.null(var_null_fit)) {
      LR_var <- df_var <- NA
    } else {
      LR_var <- as.numeric(2*logLik(alt_fit$mean) - 2*logLik(var_null_fit$mean))
      df_var <- alt_fit$disp$rank - var_null_fit$disp$rank
    }
  }

  # joint test
  if (all(exists(x = 'mean_null_formula'), exists(x = 'var_null_formula'))) {
    joint_null_fit <- tryNULL(DGLM_norm(m.form = mean_null_formula,
                                        d.form = var_null_formula,
                                        data = this_locus_df))

    if (is.null(joint_null_fit)) {
      LR_joint <- df_joint <- NA
    } else {
      LR_joint <- as.numeric(2*logLik(alt_fit$mean) - 2*logLik(joint_null_fit$mean))
      df_joint <- alt_fit$mean$rank + alt_fit$disp$rank - joint_null_fit$disp$rank - joint_null_fit$mean$rank
    }
  }

  dplyr::data_frame(LR_mean = LR_mean,
                    df_mean = df_mean,
                    LR_var = LR_var,
                    df_var = df_var,
                    LR_joint = LR_joint,
                    df_joint = df_joint)
}




fit_JRS <- function(mean_alt_formula,
                    var_alt_formula,
                    mean_null_formula,
                    var_null_formula,
                    this_locus_df) {

  return(NULL)
}