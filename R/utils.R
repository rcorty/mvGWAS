#' @title tryNA
#' @description silently get what you can out of an expression
#'
#' @param expr the expression to try to evaluate
#'
#' @return If the expression evaluates without error, the evaluated expression.
#' If there is an error in evaluating the expression, NA.
#' Warnings are suppresssed.
#'
#' @export
#'
tryNA <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NA,
                            finally = NA))
}


#' @title tryNULL
#' @description silently get what you can out of an expression
#'
#' @param expr the expression to try to evaluate
#'
#' @return If the expression evaluates without error, the evaluated expression.
#' If there is an error in evaluating the expression, NULL.
#' Warnings are suppresssed.
#'
#' @export
#'
tryNULL <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NULL,
                            finally = NULL))
}






#' @title scan_vcf_file_
#' @name mvGWAS_internals_
#'
#' @rdname mvGWAS_internals
#'
#' @return nothing
#'
#' @importFrom dplyr %>%
#'
scan_vcf_file_ <-  function(file_name,
                            drop_gts_w_fewer_than_x_obs = 5) {

  message('Started scan_vcf_file_')

  # crazy but necessary to use rslurm
  try(expr = attach(add_objects), silent = TRUE)

  vcf <- vcfR::read.vcfR(file = file_name, verbose = FALSE)

  IDs <- vcf@gt %>% colnames() %>% .[-1]
  num_snps <- dim(vcf)[1]
  num_indivs <- length(IDs)
  LR_mean <- LR_var <- LR_joint <- rep(NA, num_snps)
  df_mean <- df_var <- df_joint <- rep(NA, num_snps)
  beta_mean <- se_mean <- rep(NA, num_snps)
  beta_var <- se_var <- rep(NA, num_snps)
  this_locus_n <- rep(NA, num_snps)
  gts_raw <- rep(0, num_indivs)


  message('Started for loop over SNPs')

  for (snp_idx in 1:num_snps) {

    # drop FORMAT columns
    snp <- vcf@gt[snp_idx, -1]

    # pull out first column...may want to check in FORMAT column that MAP_gt is first
    last_gts <- gts_raw
    gts_raw <- stringr::str_split(string = snp, pattern = ':', simplify = TRUE)[,1]

    # if no genetic variation at this locus, move on
    if (length(unique(gts_raw)) == 1) { next }

    # if this snp is exactly the same as the last one, move on
    # may want to move on if it even creates the same grouping?
    if (all(last_gts == gts_raw)) {
      LR_mean[snp_idx]  <- LR_mean[snp_idx - 1]
      LR_var[snp_idx]   <- LR_var[snp_idx - 1]
      LR_joint[snp_idx] <- LR_joint[snp_idx - 1]
      df_mean[snp_idx] <- df_mean[snp_idx - 1]
      df_var[snp_idx] <- df_var[snp_idx - 1]
      df_joint[snp_idx] <- df_joint[snp_idx - 1]
      next
    }


    # if hets in both 'directions' are present, collapse them to one
    # may need to deal with '0\1' or '0/1' at some point
    if (all('0|1' %in% gts_raw, '1|0' %in% gts_raw)) {
      gts <- replace(x = gts_raw, list = gts_raw == '1|0', values = '0|1')
    } else {
      gts <- gts_raw
    }

    # if there's any very rare GT, drop it
    if (any(table(gts) < drop_gts_w_fewer_than_x_obs)) {
      bad_gts <- names(table(gts))[table(gts) < drop_gts_w_fewer_than_x_obs]
      gts <- replace(x = gts, list = which(gts == bad_gts), values = NA)
    }

    # if there's only one level, move on
    # have to check again after possibly dropping a rare GT
    if (length(unique(gts)) == 1) { next }

    # after all that checking and pruning, finally turn gts into factor
    gts <- factor(gts)

    this_locus_df <- phenotype_df
    this_locus_df[['MAP_gt']] <- gts[match(x = IDs, table = phenotype_df[[1]])]
    this_locus_df <- na.omit(object = this_locus_df)
    this_locus_n <- nrow(this_locus_df)

    alt_fit <- tryNULL(dglm::dglm(formula = mean_formula,
                                  dformula = var_formula,
                                  data = this_locus_df))

    # if we couldn't fit the alt model, move on
    if (is.null(alt_fit)) { next }

    # mean test
    if (!is.null(mean_null_formula)) {
      mean_null_fit <- tryNULL(dglm::dglm(formula = mean_null_formula,
                                          dformula = var_formula,
                                          data = this_locus_df))

      if (is.null(mean_null_fit)) { next }
      LR_mean[snp_idx] <- mean_null_fit$m2loglik - alt_fit$m2loglik
      df_mean[snp_idx] <- df.residual(mean_null_fit) - df.residual(alt_fit)
    }

    # var test
    if (!is.null(var_null_formula)) {
      var_null_fit <- tryNULL(dglm::dglm(formula = mean_formula,
                                         dformula = var_null_formula,
                                         data = this_locus_df))

      if (is.null(var_null_fit)) { next }
      LR_var[snp_idx] <- var_null_fit$m2loglik - alt_fit$m2loglik
      df_var[snp_idx] <- df.residual(var_null_fit$dispersion.fit) - df.residual(alt_fit$dispersion.fit)
    }

    # joint test
    if (!is.null(mean_null_formula) & !is.null(var_null_formula)) {
      joint_null_fit <- tryNULL(dglm::dglm(formula = mean_null_formula,
                                           dformula = var_null_formula,
                                           data = this_locus_df))

      if (is.null(joint_null_fit)) { next }
      LR_joint[snp_idx] <- joint_null_fit$m2loglik - alt_fit$m2loglik
      df_joint[snp_idx] <- df.residual(mean_null_fit) - df.residual(alt_fit) +
        df.residual(var_null_fit$dispersion.fit) - df.residual(alt_fit$dispersion.fit)
    }

  }
  message('Finished for loop over SNPs')

  fix_df <- vcf@fix %>%
    dplyr::as_data_frame() %>%
    dplyr::select(-QUAL, -INFO) %>%
    dplyr::mutate(POS = as.numeric(POS))

  result <- dplyr::data_frame(LR_mean = LR_mean,
                              LR_var = LR_var,
                              LR_joint = LR_joint,
                              df_mean = df_mean,
                              df_var = df_var,
                              df_joint = df_joint,
                              mean_asymp_p  = pchisq(q = LR_mean,  df = df_mean,  lower.tail = FALSE),
                              var_asymp_p   = pchisq(q = LR_var,   df = df_var,   lower.tail = FALSE),
                              joint_asymp_p = pchisq(q = LR_joint, df = df_joint, lower.tail = FALSE))


  return(dplyr::bind_cols(fix_df, result))
}

