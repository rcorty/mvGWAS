#' @title mvGWAS reference class
#' @name mvGWAS
#'
#' @description A reference class implementing a Mean-Variance Genome Wide Association Study
#'
#' @field metadata information about the mvGWAS
#' @field data the data used to conduct the mvGWAS
#' @field intermediates computations stored in the course of conducting the mvGWAS that aren't readily interpretable
#' @field results the results of the mvGWAS
#'
#' @export mvGWAS
#' @exportClass mvGWAS
#'
mvGWAS <- setRefClass(
  Class = 'mvGWAS',
  fields = list(metadata = 'list',
                data = 'list',
                intermediataes = 'list',
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
    if (!(names(phenotypes)[1] %in% c('ID', 'IID'))) {
      warning("First column of phenotype file will be used as ID column for matching genotypes.
                 It's not named 'ID' nor 'IID', so make sure it's there!")
    }

    # weakly validate genotype directory
    # definitely a dir containing at least ont .VCF.GZ file
    # maybe it can contain only that type of files?


    metadata <<- list(created_at = Sys.time(),
                      phenotype_file = phenotype_file,
                      genotype_directory = genotype_directory)

    data <<- list(phenotypes = phenotypes)

  }
)



#' @title scan_vcf_file_
#' @name mvGWAS_internals_
#'
#' @rdname mvGWAS_internals
#'
#' @return nothing
#'
mvGWAS$methods(
  scan_vcf_file_ = function(file_name,
                            drop_gts_w_fewer_than_x_obs = 5) {

    vcf <- vcfR::read.vcfR(file = file_name, verbose = FALSE)

    IDs <- vcf@gt %>% colnames() %>% .[-1]
    num_snps <- dim(vcf)[1]
    num_indivs <- length(IDs)
    LR_mean <- LR_var <- LR_joint <- rep(NA, num_snps)
    df_mean <- df_var <- df_joint <- rep(NA, num_snps)
    gts_raw <- rep(0, num_indivs)

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

      this_locus_df <- data$phenotypes
      this_locus_df[['MAP_gt']] <- gts[match(x = IDs, table = data$phenotypes[[1]])]
      this_locus_df <- na.omit(object = this_locus_df)

      alt_fit <- tryNULL(dglm::dglm(formula = metadata$mean_formula,
                                    dformula = metadata$var_formula,
                                    data = this_locus_df))

      # if we couldn't fit the alt model, move on
      if (is.null(alt_fit)) { next }

      # mean test
      if ('mean_null_formula' %in% names(metadata)) {
        mean_null_fit <- tryNULL(dglm::dglm(formula = metadata$mean_null_formula,
                                            dformula = metadata$var_formula,
                                            data = this_locus_df))

        if (is.null(mean_null_fit)) { next }
        LR_mean[snp_idx] <- mean_null_fit$m2loglik - alt_fit$m2loglik
        df_mean[snp_idx] <- df.residual(mean_null_fit) - df.residual(alt_fit)
      }

      # var test
      if ('var_null_formula' %in% names(metadata)) {
        var_null_fit <- tryNULL(dglm::dglm(formula = metadata$mean_formula,
                                           dformula = metadata$var_null_formula,
                                           data = this_locus_df))

        if (is.null(var_null_fit)) { next }
        LR_var[snp_idx] <- var_null_fit$m2loglik - alt_fit$m2loglik
        df_var[snp_idx] <- df.residual(var_null_fit$dispersion.fit) - df.residual(alt_fit$dispersion.fit)
      }

      # joint test
      if (all(c('mean_null_formula', 'var_null_formula') %in% names(metadata))) {
        joint_null_fit <- tryNULL(dglm::dglm(formula = metadata$mean_null_formula,
                                             dformula = metadata$var_null_formula,
                                             data = this_locus_df))

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
    stopifnot(all.vars(expr = formula) %in% c(names(data$phenotypes), 'MAP_gt'))  # todo: add other gt info -- dosage and genoprobs

    metadata[[paste0(type, '_formula')]] <<- formula

    rhs <- formula[[switch(EXPR = type, mean = 3, var = 2)]]
    rhs_vars <- all.vars(rhs)
    special_term_idxs <- grep(pattern = 'MAP_gt', x = rhs_vars)

    if (any(special_term_idxs)) {

      # corner case where there are no non-special RHS terms
      if (length(special_term_idxs) == length(rhs_vars)) {

        null_formula <- stats::reformulate(termlabels = '1',
                                           response = switch(EXPR = type,
                                                             mean = formula[[2]],
                                                             var = NULL))

        # typical case
      } else {

        null_formula <- stats::formula(x = stats::drop.terms(termobj = stats::terms(formula),
                                                             dropx = special_term_idxs,
                                                             keep.response = switch(EXPR = type,
                                                                                    mean = TRUE,
                                                                                    var = FALSE)))

      }

      metadata[[paste0(type, '_null_formula')]] <<- null_formula

    }
  }
)



#' @title conduct_scan_local_
#' @name mvGWAS_internals
#'
#' @param num_cores the number of cores to use, defaults to parallel::detectCores()
#'
#' @description conducts a genome scan locally
#'
#' @return nothing
#'
mvGWAS$methods(
  conduct_scan_local_ = function(num_cores = parallel::detectCores()) {

    usingMethods(scan_vcf_file_)

    genotype_files <- list.files(path = metadata$genotype_directory, full.names = TRUE, pattern = '\\.vcf|\\.VCF')

    #this should be mclapply
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
#' @name mvGWAS_internals
#'
#' @description conducts a genome scan on a slurm cluster
#'
#' @return nothing
#'
mvGWAS$methods(
  conduct_scan_slurm_ = function(max_num_nodes = Inf) {

    usingMethods(scan_vcf_file_)

    genotype_files <- list.files(path = metadata$genotype_directory, full.names = TRUE, pattern = '\\.vcf|\\.VCF')

    sjob <- rslurm::slurm_apply(f = scan_vcf_file_,
                                params = dplyr::data_frame(file_name = genotype_files),
                                nodes = min(max_num_nodes, length(genotype_files)))

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

    usingMethods(stash_formula_, conduct_scan_local_, conduct_scan_slurm_)

    system <- match.arg(arg = system)

    stash_formula_(formula = mean_formula, type = 'mean')
    stash_formula_(formula = var_formula, type = 'var')

    if (!any(grepl(pattern = 'null_formula', x = names(metadata)))) {
      stop('Use MAP_gt in at least one of mean_formula or var_formula.')
    }

    switch(EXPR = system,
           'local' = conduct_scan_local_(...),
           'slurm' = conduct_scan_slurm_(...))

    return(.self)

  }
)




