#' @title A Reference Class to represent a mean-variance GWAS
#' @name mvGWAS
#'
#' @description A reference class implementing a Mean-Variance Genome Wide Association Study
#'
#' @field metadata information about the mvGWAS
#' @field data the data used to conduct the mvGWAS
#' @field null_model the model fit without any SNPs
#' @field results the results of the mvGWAS
#'
#' @export mvGWAS
#' @exportClass mvGWAS
#'
mvGWAS <- setRefClass(
  Class = 'mvGWAS',
  fields = list(metadata = 'list',
                data = 'list',
                null_model = 'ANY',    # todo: make this DGLM, but need to import it beforehand or something
                results = 'data.frame'))


#' @title initialize
#' @name mvGWAS_initialize
#'
#' @description Initialize an mvGWAS object
#'
#' @param phenotype_file The file that stores the phenotypes.  Must be a .CSV or .RDS.  Must have one column called 'ID'.
#' @param genotype_directory The directory that stores the genotypes.  They must be in .vcf.gz format.
#'
#' @return an mvGWAS object
#'
#' @importFrom dplyr %>%
#'
mvGWAS$methods(
  initialize = function(phenotype_file,
                        genotype_directory) {

    # check inputs
    stopifnot(file.exists(phenotype_file))
    stopifnot(dir.exists(genotype_directory))

    # read in phenotypes
    if (phenotype_file %>% tools::file_ext() %>% toupper() == 'RDS') {

      phenotypes <- readRDS(file = phenotype_file)

    } else if (phenotype_file %>% tools::file_ext() %>% toupper() == 'CSV') {

      phenotypes <- readr::read_csv(phenotype_file)

    } else {

      stop('phenotype_file must be a .RDS or .CSV')

    }

    # a little checking of phenotype file
    if (names(phenotypes)[1] != 'ID') {
      stop('First column of phenotype_file must be "ID".  It will be used for matching with genotypes.')
    }
    if (nrow(phenotypes) > 1e5) {
      warning('mvGWAS hasnt been tested on datasets with more than 100,000 individuals.')
    }

    # weakly validate genotype directory
    # definitely a dir containing at least ont .VCF.GZ file
    # maybe it can contain only that type of files?


    metadata <<- list(created_at = Sys.time(),
                      phenotype_file = phenotype_file,
                      genotype_directory = genotype_directory,
                      available_keywords = c('GT', 'DS', 'GP'))

    data <<- list(phenotypes = phenotypes)

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
        new_terms <- gsub(pattern = 'GP', replacement = 'GP_het + GP_alt', x = labels(stats::terms(formula)))
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




#' @title scan_vcf_file_
#' @name mvGWAS_internals
#'
#' @return nothing
#'
#' @importFrom dplyr %>%
#'
mvGWAS$methods(
  scan_vcf_file_ = function(file_name,
                            min_gt_count = 5,
                            min_gp = 0.05) {

    try(expr = attach(what = add_objects), silent = TRUE)
    try(expr = attach(what = data, warn.conflicts = FALSE), silent = TRUE)
    try(expr = attach(what = metadata, warn.conflicts = FALSE), silent = TRUE)

    vcf <- vcfR::read.vcfR(file = file_name, verbose = FALSE)

    num_snps <- dim(vcf)[1]
    num_indivs <- dim(vcf)[3] - 1

    LR_mean <- LR_var <- LR_joint <- rep(NA, num_snps)
    df_mean <- df_var <- df_joint <- rep(NA, num_snps)
    this_locus_n <- rep(NA, num_snps)

    if ('DS' %in% all.vars(mean_alt_formula)) {
      DS_beta_mean <- DS_se_mean <- rep(NA, num_snps)
    }
    if ('DS' %in% all.vars(var_alt_formula)) {
      DS_beta_var <- DS_se_var <- rep(NA, num_snps)
    }

    last_locus_df <- dplyr::data_frame()

    for (snp_idx in 1:num_snps) {

      this_locus_df <- phenotypes

      # drop FORMAT columns
      snp_row <- vcf@gt[snp_idx, ]

      if ('GT' %in% used_keywords) {
        this_locus_df <- dplyr::inner_join(x = this_locus_df,
                                           y = pull_GT(snp_row = snp_row, min_gt_count = min_gt_count),
                                           by = 'ID')
      }
      if ('DS' %in% used_keywords) {
        this_locus_df <- dplyr::inner_join(x = this_locus_df,
                                           y = pull_DS(snp_row = snp_row),
                                           by = 'ID')
      }
      if ('GP' %in% used_keywords) {
        this_locus_df <- dplyr::inner_join(x = this_locus_df,
                                           y = pull_GP(snp_row = snp_row, min_gp = min_gp),
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

      alt_fit <- tryNULL(dglm::dglm(formula = mean_alt_formula,
                                    dformula = var_alt_formula,
                                    data = this_locus_df,
                                    method = 'ml'))

      # if we couldn't fit the alt model, move on
      if (is.null(alt_fit)) { next }

      if ('DS' %in% all.vars(mean_alt_formula)) {
        DS_beta_mean[snp_idx] <- coef(summary(alt_fit))['DS', 'Estimate']
        DS_se_mean[snp_idx] <- coef(summary(alt_fit))['DS', 'Std. Error']
      }
      if ('DS' %in% all.vars(var_alt_formula)) {
        DS_beta_var[snp_idx] <- coef(summary(alt_fit$dispersion.fit))['DS', 'Estimate']
        DS_se_var[snp_idx] <-coef(summary(alt_fit$dispersion.fit))['DS', 'Std. Error']
      }

      # mean test
      if (exists(x = 'mean_null_formula')) {
        mean_null_fit <- tryNULL(dglm::dglm(formula = mean_null_formula,
                                            dformula = var_alt_formula,
                                            data = this_locus_df,
                                            method = 'ml'))

        if (is.null(mean_null_fit)) { next }
        LR_mean[snp_idx] <- mean_null_fit$m2loglik - alt_fit$m2loglik
        df_mean[snp_idx] <- df.residual(mean_null_fit) - df.residual(alt_fit)
      }

      # var test
      if (exists(x = 'var_null_formula')) {
        var_null_fit <- tryNULL(dglm::dglm(formula = mean_alt_formula,
                                           dformula = var_null_formula,
                                           data = this_locus_df,
                                           method = 'ml'))

        if (is.null(var_null_fit)) { next }
        LR_var[snp_idx] <- var_null_fit$m2loglik - alt_fit$m2loglik
        df_var[snp_idx] <- df.residual(var_null_fit$dispersion.fit) - df.residual(alt_fit$dispersion.fit)
      }

      # joint test
      if (all(exists(x = 'mean_null_formula'), exists(x = 'var_null_formula'))) {
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
      dplyr::mutate(POS = as.integer(POS))

    result <- dplyr::data_frame(n             = this_locus_n,
                                LR_mean       = LR_mean,
                                LR_var        = LR_var,
                                LR_joint      = LR_joint,
                                df_mean       = df_mean,
                                df_var        = df_var,
                                df_joint      = df_joint,
                                mean_asymp_p  = pchisq(q = LR_mean,  df = df_mean,  lower.tail = FALSE),
                                var_asymp_p   = pchisq(q = LR_var,   df = df_var,   lower.tail = FALSE),
                                joint_asymp_p = pchisq(q = LR_joint, df = df_joint, lower.tail = FALSE))


    if ('DS' %in% all.vars(mean_alt_formula)) {
      result <- dplyr::bind_cols(result, DS_beta_mean = DS_beta_mean, DS_se_mean = DS_se_mean, DS_z_mean = DS_beta_mean/DS_se_mean)
    }
    if ('DS' %in% all.vars(var_alt_formula)) {
      result <- dplyr::bind_cols(result, DS_beta_var = DS_beta_var, DS_se_var = DS_se_var, DS_z_var = DS_beta_var/DS_se_var)
    }

    return(dplyr::bind_cols(fix_df, result))
  }
)


#' @title conduct_scan_local_
#' @name mvGWAS_conduct_scan
#'
#' @param num_cores the number of cores to use, defaults to parallel::detectCores()
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

    genotype_files <- list.files(path = metadata$genotype_directory, full.names = TRUE, pattern = '\\.vcf|\\.VCF')

    # better to use lapply directly if num_cores == 1 rather than let parallel::mclapply call it
    # bc user likely doesn't have parallel package installed
    if (num_cores == 1) {
      results_list <- lapply(X = genotype_files,
                             FUN = scan_vcf_file_)
    } else {
      results_list <- parallel::mclapply(X = genotype_files,
                                         FUN = scan_vcf_file_)
    }

    results <<- results_list %>% dplyr::bind_rows()
  }
)


#' @title conduct_scan_slurm_
#' @name mvGWAS_conduct_scan
#'
#' @description conducts a genome scan on a slurm cluster
#'
#' @return nothing
#'
#' @importFrom dplyr %>%
#'
mvGWAS$methods(
  conduct_scan_slurm_ = function(max_num_nodes = Inf, ...) {

    usingMethods(scan_vcf_file_)

    genotype_files <- list.files(path = metadata$genotype_directory, full.names = TRUE, pattern = '\\.vcf|\\.VCF')

    sjob <- rslurm::slurm_apply(f = scan_vcf_file_,
                                params = dplyr::data_frame(file_name = genotype_files),
                                nodes = min(max_num_nodes, length(genotype_files)),
                                cpus_per_node = 1,
                                add_objects = list(mean_alt_formula = .self$metadata$mean_alt_formula,
                                                   var_alt_formula = .self$metadata$var_alt_formula,
                                                   mean_null_formula = .self$metadata$mean_null_formula,
                                                   var_null_formula = .self$metadata$var_null_formula,
                                                   used_keywords = .self$metadata$used_keywords,
                                                   phenotypes = .self$data$phenotypes))


    results_list <- rslurm::get_slurm_out(slr_job = sjob, outtype = "raw", wait = TRUE)

    results <<- results_list %>% dplyr::bind_rows()

  }
)


#' @title conduct_scan
#' @name mvGWAS_conduct_scan
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

    usingMethods(determine_keyword_use_, stash_formula_, conduct_scan_local_, conduct_scan_slurm_)

    system <- match.arg(arg = system)

    determine_keyword_use_(mean_formula, var_formula)

    stash_formula_(formula = mean_formula, type = 'mean')
    stash_formula_(formula = var_formula, type = 'var')

    null_model <<- dglm::dglm(formula = metadata$mean_null_formula,
                              dformula = metadata$var_null_formula,
                              data = data$phenotypes)

    if (length(metadata$used_keywords) == 0) {
      message('No keywords in mean_formula nor var_formula.  Returning without scanny any SNPs.')
      return(TRUE)
    }

    switch(EXPR = system,
           'local' = conduct_scan_local_(...),
           'slurm' = conduct_scan_slurm_(...))

    return(TRUE)

  }
)




#' @title write_results_to_file
#' @name mvGWAS_write_results_to_file
#'
#' @description Writes results to file
#'
#' @return TRUE if successful, FALSE otherwise.
#' @importFrom dplyr %>%
#'
mvGWAS$methods(
  write_results_to_file = function(file_name) {

    if (missing(file_name)) {
      phenotype_name <- as.character(x = metadata$mean_alt_formula[[2]])
      file_name <- paste0(phenotype_name, '_', gsub(pattern = ' ', replacement = '_', x = Sys.time()), '.tsv')
    }

    readr::write_tsv(x = results %>% mutate_if(is.double, round, 3),
                     path = file_name,
                     col_names = TRUE)

  }
)






#' @title manhattan_plot
#' @name mvGWAS_manhattan_plot
#'
#' @description makes a Manhattan plot
#'
#' @return the plot
#' @importFrom dplyr %>%
#'
mvGWAS$methods(
  manhattan_plot = function(what_to_plot = c('asymp_p', 'asymp_p_gc'),
                            lowest_nlog_p_to_plot = 3,
                            highlight_nlog_p_above = 5,
                            highlight_size = 2,
                            nonhighlight_alpha = 0.5) {

    what_to_plot <- match.arg(arg = what_to_plot)

    if (what_to_plot == 'asymp_p') {

      for_plotting <- results %>%
        tidyr::gather(key = test, value = p, mean_asymp_p, var_asymp_p, joint_asymp_p) %>%
        dplyr::mutate(test = dplyr::recode(test, mean_asymp_p = 'mean', var_asymp_p = 'var', joint_asymp_p = 'joint'))

    } else {

      stop('mahnattan_plot() on gc values not yet implemented')

    }

    for_plotting2 <- for_plotting %>%
      dplyr::mutate(nlog_p = -log(p, base = 10)) %>%
      dplyr::filter(nlog_p > lowest_nlog_p_to_plot, nlog_p < highlight_nlog_p_above)

    for_plotting3 <- for_plotting %>%
      dplyr::mutate(nlog_p = -log(p, base = 10)) %>%
      dplyr::filter(nlog_p > highlight_nlog_p_above)

    ggplot2::ggplot(data = for_plotting2, mapping = ggplot2::aes(x = POS, y = nlog_p, color = test)) +
      ggplot2::facet_grid(~CHROM, scales = 'free_x', switch = 'x', space = 'free', drop = FALSE) +
      ggplot2::geom_point(alpha = nonhighlight_alpha) +
      ggplot2::geom_point(data = for_plotting3, size = highlight_size) +
      ggplot2::scale_color_manual(values = c(mean = 'blue', var = 'red', joint = 'black')) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank()) +
      ggplot2::xlab('Chromosome and position in Mb') +
      ggplot2::ylab('-log10(p)')

  }
)




#' @title qq_plot
#' @name mvGWAS_qq_plot
#'
#' @description makes a Manhattan plot
#'
#' @return the plot
#' @importFrom dplyr %>%
#'
mvGWAS$methods(
  qq_plot = function(what_to_plot,
                     theoretical = 'unif') {

    v <- results[[what_to_plot]]

    for_plotting <- dplyr::data_frame(obs = sort(v),
                                      theo = seq(from = 0, to = 1, length.out = sum(!is.na(v))))

    ggplot2::ggplot() +
      ggplot2::geom_segment(data = dplyr::data_frame(start = 0, stop = 1),
                            mapping = ggplot2::aes(x = start, y = start, xend = stop, yend = stop),
                            color = 'gray') +
      ggplot2::geom_point(data = for_plotting,
                          mapping = ggplot2::aes(x = theo, y = obs)) +
      ggplot2::ggtitle(label = paste0('QQ plot of ', what_to_plot)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank()) +
      ggplot2::xlab('Theoretical uniform distribution') +
      ggplot2::ylab('Observed p values')


  }
)


