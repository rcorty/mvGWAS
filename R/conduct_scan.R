#' @title scan_vcf_file_
#' @name mvGWAS_internals
#'
#' @return nothing
#' @importFrom dplyr %>%
#'
mvGWAS$methods(

  scan_vcf_file_ = function(file_name,
                            min_gt_count = 5) {

    # for local and slurm
    # try(expr = attach(what = data, warn.conflicts = FALSE), silent = TRUE)
    # try(expr = attach(what = metadata, warn.conflicts = FALSE), silent = TRUE)

    # for SLURM
    # try(expr = attach(what = add_objects, warn.conflicts = FALSE), silent = TRUE)


    vcf <- vcfR::read.vcfR(file = file_name, verbose = FALSE)

    num_snps <- dim(vcf)[1]
    num_indivs <- dim(vcf)[3] - 1
    message('num_snps: ', num_snps)
    message('num_indivs: ', num_indivs)

    LR_mean <- LR_var <- LR_joint <- rep(NA, num_snps)
    df_mean <- df_var <- df_joint <- rep(NA, num_snps)
    this_locus_n <- rep(NA, num_snps)

    if ('DS' %in% all.vars(metadata$mean_alt_formula)) {
      beta_DS_mean <- se_DS_mean <- rep(NA, num_snps)
    }
    if ('DS' %in% all.vars(metadata$var_alt_formula)) {
      beta_DS_var <- se_DS_var <- rep(NA, num_snps)
    }

    last_locus_df <- dplyr::data_frame()

    for (snp_idx in 1:num_snps) {

      # message('Starting SNP ', snp_idx)

      this_locus_df <- data$phenotypes

      # drop FORMAT columns
      snp_row <- vcf@gt[snp_idx, ]

      if ('GT' %in% metadata$used_keywords) {
        this_locus_df <- dplyr::inner_join(x = this_locus_df,
                                           y = pull_GT(snp_row = snp_row, min_gt_count = min_gt_count),
                                           by = 'ID')
      }
      if ('DS' %in% metadata$used_keywords) {
        this_locus_df <- dplyr::inner_join(x = this_locus_df,
                                           y = pull_DS(snp_row = snp_row),
                                           by = 'ID')
      }
      if ('GP' %in% metadata$used_keywords) {
        this_locus_df <- dplyr::inner_join(x = this_locus_df,
                                           y = pull_GP(snp_row = snp_row),
                                           by = 'ID')
      }

      this_locus_df <- na.omit(object = this_locus_df)

      if (identical(this_locus_df, last_locus_df)) {
        LR_mean[snp_idx]  <- LR_mean[snp_idx - 1]
        df_mean[snp_idx]  <- df_mean[snp_idx - 1]
        LR_var[snp_idx]   <- LR_var[snp_idx - 1]
        df_var[snp_idx]   <- df_var[snp_idx - 1]
        LR_joint[snp_idx] <- LR_joint[snp_idx - 1]
        df_joint[snp_idx] <- df_joint[snp_idx - 1]
        next
      }

      this_locus_n[snp_idx] <- nrow(this_locus_df)
      if (nrow(this_locus_df) == 0) { next }

      alt_fit <- tryNULL(dglm::dglm(formula = metadata$mean_alt_formula,
                                    dformula = metadata$var_alt_formula,
                                    data = this_locus_df,
                                    method = 'ml'))

      # if we couldn't fit the alt model, COMPLETELY, move on
      # I believe that when one coef is NA (not estimated) it causes problems later
      # e.g. with calculating the df of a test
      if (is.null(alt_fit)) {
        next
      }
      if (any(is.na(coef(alt_fit)), is.na(coef(alt_fit$dispersion.fit)))) {
        next
      }

      if ('DS' %in% all.vars(metadata$mean_alt_formula)) {

        beta_DS_mean[snp_idx] <- tryNA(coef(summary(alt_fit))['DS', 'Estimate'])
        se_DS_mean[snp_idx] <- tryNA(coef(summary(alt_fit))['DS', 'Std. Error'])
      }

      if ('DS' %in% all.vars(metadata$var_alt_formula)) {

        beta_DS_var[snp_idx] <- tryNA(coef(summary(alt_fit$dispersion.fit))['DS', 'Estimate'])
        se_DS_var[snp_idx] <- tryNA(coef(summary(alt_fit$dispersion.fit))['DS', 'Std. Error'])
      }

      # mean test
      if (exists(x = 'mean_null_formula', where = metadata)) {
        mean_null_fit <- tryNULL(dglm::dglm(formula = mean_null_formula,
                                            dformula = var_alt_formula,
                                            data = this_locus_df,
                                            method = 'ml'))

        if (is.null(mean_null_fit)) { next }
        LR_mean[snp_idx] <- mean_null_fit$m2loglik - alt_fit$m2loglik
        df_mean[snp_idx] <- df.residual(mean_null_fit) - df.residual(alt_fit)

        # if (df_mean[snp_idx] == 0) { browser() }
      }

      # var test
      if (exists(x = 'var_null_formula', where = metadata)) {
        var_null_fit <- tryNULL(dglm::dglm(formula = mean_alt_formula,
                                           dformula = var_null_formula,
                                           data = this_locus_df,
                                           method = 'ml'))

        if (is.null(var_null_fit)) { next }
        LR_var[snp_idx] <- var_null_fit$m2loglik - alt_fit$m2loglik
        df_var[snp_idx] <- df.residual(var_null_fit$dispersion.fit) - df.residual(alt_fit$dispersion.fit)
      }

      # joint test
      if (all(exists(x = 'metadata$mean_null_formula', where = metadata),
              exists(x = 'var_null_formula', where = metadata))) {
        joint_null_fit <- tryNULL(dglm::dglm(formula = mean_null_formula,
                                             dformula = var_null_formula,
                                             data = this_locus_df,
                                             method = 'ml'))

        if (is.null(joint_null_fit)) { next }
        LR_joint[snp_idx] <- joint_null_fit$m2loglik - alt_fit$m2loglik
        df_joint[snp_idx] <- df.residual(mean_null_fit) +
          df.residual(var_null_fit$dispersion.fit) -
          df.residual(alt_fit) -
          df.residual(alt_fit$dispersion.fit)
      }

    }

    fix_df <- vcf@fix %>%
      dplyr::as_data_frame() %>%
      dplyr::select(-QUAL, -INFO) %>%
      dplyr::mutate(POS = as.integer(POS),
                    vcf_file = file_name)

    result <- dplyr::data_frame(n          = this_locus_n,
                                LR_mean    = LR_mean,
                                LR_var     = LR_var,
                                LR_joint   = LR_joint,
                                df_mean    = df_mean,
                                df_var     = df_var,
                                df_joint   = df_joint,
                                p_LR_mean  = pchisq(q = LR_mean,  df = df_mean,  lower.tail = FALSE),
                                p_LR_var   = pchisq(q = LR_var,   df = df_var,   lower.tail = FALSE),
                                p_LR_joint = pchisq(q = LR_joint, df = df_joint, lower.tail = FALSE))

    if ('DS' %in% all.vars(metadata$mean_alt_formula)) {
      result <- dplyr::bind_cols(result,
                                 dplyr::data_frame(beta_DS_mean = beta_DS_mean,
                                                   se_DS_mean = se_DS_mean,
                                                   z_DS_mean = beta_DS_mean/se_DS_mean,
                                                   p_z_DS_mean = dplyr::if_else(condition = z_DS_mean < 0,
                                                                                true = pnorm(q = z_DS_mean),
                                                                                false = pnorm(q = z_DS_mean, lower.tail = TRUE),
                                                                                missing = NA_real_)))
    }
    if ('DS' %in% all.vars(metadata$var_alt_formula)) {
      result <- dplyr::bind_cols(result,
                                 dplyr::data_frame(beta_DS_var = beta_DS_var,
                                                   se_DS_var = se_DS_var,
                                                   z_DS_var = beta_DS_var/se_DS_var,
                                                   p_z_DS_var = dplyr::if_else(condition = z_DS_var < 0,
                                                                               true = pnorm(q = z_DS_var),
                                                                               false = pnorm(q = z_DS_var, lower.tail = TRUE),
                                                                               missing = NA_real_)))
    }

    return(dplyr::bind_cols(fix_df, result))
  }

)



#' @title determine_keyword_use_
#' @name mvGWAS_internals
#'
#' @description determines which keywords are used in mean_formula and var_formula
#'
#' @return nothing
#'
mvGWAS$methods(

  determine_keyword_use_ = function(mean_formula, var_formula) {

    used_keyword_idxs <- metadata$available_keywords %in% c(all.vars(mean_formula), all.vars(var_formula))

    metadata$used_keywords <<- metadata$available_keywords[used_keyword_idxs]
  }

)


#' @title stash_formula_
#' @name mvGWAS_internals
#'
#' @param formula the formula to validate
#' @param type one of 'mean' or 'var'
#'
#' @description validates a formula and stashes it for use in an mvGWAS
#'
#' @return nothing
#'
mvGWAS$methods(

  stash_formula_ = function(formula, type = c('mean', 'var')) {

    type <- match.arg(arg = type)

    stopifnot(inherits(x = formula, 'formula'))
    stopifnot(length(formula) == switch(EXPR = type, mean = 3, var = 2))
    stopifnot(all.vars(expr = formula) %in% c(names(data$phenotypes), metadata$available_keywords))

    # derive and stash null formula
    rhs <- formula[[switch(EXPR = type, mean = 3, var = 2)]]
    rhs_vars <- all.vars(rhs)
    keyword_idxs <- grep(pattern = paste(metadata$available_keywords, collapse = '|'), x = rhs_vars)

    # If the formula uses one or more keywords (the words used to refer to the locus),
    # remove them to create a null formula then save the original as the alternative formula
    if (any(keyword_idxs)) {

      # corner case where there are no non-keyword RHS terms
      if (length(keyword_idxs) == length(rhs_vars)) {

        r <- switch(EXPR = type, mean = formula[[2]], var = NULL)
        null_formula <- stats::reformulate(termlabels = '1', response = r)

        # typical case
      } else {

        kr <- switch(EXPR = type, mean = TRUE, var = FALSE)
        keyword_free <- stats::drop.terms(termobj = stats::terms(formula),
                                          dropx = keyword_idxs,
                                          keep.response = kr)
        null_formula <- stats::formula(x = keyword_free)

      }

      # modify formula if it uses GP -- need two terms -- one for het and one for homozygous alternative
      if ('GP' %in% all.vars(expr = formula)) {

        r <- switch(EXPR = type, mean = formula[[2]], var = NULL)
        new_terms <- gsub(pattern = 'GP', replacement = 'GP_add + GP_dom', x = labels(stats::terms(formula)))
        alt_formula <- reformulate(termlabels = new_terms, response = r)

      } else {
        alt_formula <- formula
      }

      # If the formula doesn't use any keywords, just save it as the null and we don't have an alt
    } else {

      alt_formula <- NULL
      null_formula <- formula

    }

    metadata[[paste0(type, '_null_formula')]] <<- null_formula
    metadata[[paste0(type, '_alt_formula')]] <<- alt_formula
  }

)





#' @title conduct_scan
#' @name mvGWAS_conduct_scan
#'
#' @param num_cores the number of cores to use for a local scan.  Defaults to \code{parallel::detectCores()}.
#'
#' @description conducts a genome scan locally
#'
#' @return nothing
#'
#' @importFrom dplyr %>%
#'
mvGWAS$methods(

  conduct_scan_local_ = function(num_cores = parallel::detectCores()) {

    usingMethods(scan_vcf_file_)

    genotype_files <- list.files(path = metadata$genotype_directory,
                                 full.names = TRUE,
                                 pattern = metadata$genotype_file_pattern)

    # better to use lapply directly if num_cores == 1 rather than let parallel::mclapply call it
    # bc user likely doesn't have parallel package installed
    if (num_cores == 1) {
      results_list <- lapply(X = genotype_files,
                             FUN = scan_vcf_file_)
    } else {
      results_list <- parallel::mclapply(X = genotype_files,
                                         FUN = scan_vcf_file_)
    }

    results <<- results_list %>%
      dplyr::bind_rows()
  }

)


#' @title conduct_scan
#' @name mvGWAS_conduct_scan
#'
#' @description conducts a genome scan on a slurm cluster
#'
#' @param job_time_in_mins amount of time to request from SLURM for each job (in minutes).  Defaults to 360.
#' @param vcf_files_per_job number of vcf files each job should scan.  Defaults to 1.
#' @param max_num_jobs the maximum number of jobs to submit to SLURM.
#'
#' @details The actual number of jobs submitted to the cluster is \code{min(max_num_jobs, ceiling(num_vcf_files/vcf_files_per_job))}.
#'
#' @return nothing
#'
#' @importFrom dplyr %>%
#'
mvGWAS$methods(

  conduct_scan_slurm_ = function(job_time_in_mins = 360,
                                 vcf_files_per_job = 1,
                                 max_num_jobs = Inf) {

    usingMethods(scan_vcf_file_)

    genotype_files <- list.files(path = metadata$genotype_directory,
                                 full.names = TRUE,
                                 pattern = metadata$genotype_file_pattern)

    sjob <- rslurm::slurm_apply(f = scan_vcf_file_,
                                params = dplyr::data_frame(file_name = genotype_files),
                                nodes = min(max_num_jobs, ceiling(length(genotype_files)/vcf_files_per_job)),
                                cpus_per_node = 1,
                                slurm_options = list(time = job_time_in_mins))

    results_list <- rslurm::get_slurm_out(slr_job = sjob, outtype = "raw", wait = TRUE)

    results <<- results_list %>%
      dplyr::bind_rows()
  }

)



#' @title conduct_scan
#' @name mvGWAS_conduct_scan
#' @include mvGWAS.R
#'
#' @param mean_formula the mean formula
#' @param var_formula the variance formula
#' @param computer Either 'local' (default) or 'slurm'. Should the scan be conducted on the local machine or on a slurm cluster.
#'
#' @description conducts an mvGWAS
#'
#' @return the mvGWAS object  (change this?)
#'
mvGWAS$methods(

  conduct_scan = function(mean_formula,
                          var_formula,
                          system = c('local', 'slurm'),
                          ...) {

    system <- match.arg(arg = system)

    determine_keyword_use_(mean_formula, var_formula)

    stash_formula_(formula = mean_formula, type = 'mean')
    stash_formula_(formula = var_formula, type = 'var')

    na_free_data <- na.omit(data$phenotypes[,unique(c(all.vars(metadata$mean_null_formula), all.vars(metadata$var_null_formula)))])

    null_model <<- dglm::dglm(formula = metadata$mean_null_formula,
                              dformula = metadata$var_null_formul,
                              data = na_free_data,
                              method = 'ml')

    if (length(metadata$used_keywords) == 0) {
      message('No keywords in mean_formula nor var_formula.  Returning without scanning any any SNPs.')
      return(TRUE)
    }

    switch(EXPR = system,
           'local' = conduct_scan_local_(...),
           'slurm' = conduct_scan_slurm_(...))

    results <<- results %>%
      dplyr::mutate(CHROM = factor(x = stringr::str_pad(string = CHROM, width = 2, pad = '0'),
                                   levels = sort(stringr::str_pad(string = unique(CHROM), width = 2, pad = '0'))))
    return(TRUE)
  }

)







