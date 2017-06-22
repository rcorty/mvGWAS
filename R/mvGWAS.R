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
    if (!(names(phenotypes)[1] %in% c('ID', 'IID'))) {
      warning("First column of phenotype file will be used as ID column for matching genotypes.
                 It's not named 'ID' nor 'IID', so make sure it's there!")
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

    metadata[[paste0(type, '_formula')]] <<- formula

    rhs <- formula[[switch(EXPR = type, mean = 3, var = 2)]]
    rhs_vars <- all.vars(rhs)
    keyword_idxs <- grep(pattern = paste(metadata$available_keywords, collapse = '|'), x = rhs_vars)

    # remove keywords (the words used to refer to the locus) to create null formula
    if (any(keyword_idxs)) {

      # corner case where there are no non-special RHS terms
      if (length(keyword_idxs) == length(rhs_vars)) {

        r <- switch(EXPR = type, mean = formula[[2]], var = NULL)
        null_formula <- stats::reformulate(termlabels = '1', response = r)

      # typical case
      } else {

        kr <- switch(EXPR = type, mean = TRUE, var = FALSE)
        null_formula <- stats::formula(x = stats::drop.terms(termobj = stats::terms(formula),
                                                             dropx = keyword_idxs,
                                                             keep.response = kr))

      }
    } else {

      null_formula <- NULL

    }

    metadata[[paste0(type, '_null_formula')]] <<- null_formula

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
  determine_keyword_use_ = function() {

    used_keyword_idxs <- metadata$available_keywords %in% c(all.vars(metadata$mean_formula), all.vars(metadata$var_formula))
    metadata$used_keywords <<- metadata$available_keywords[used_keyword_idxs]

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
#' @name mvGWAS_internals
#'
#' @description conducts a genome scan on a slurm cluster
#'
#' @return nothing
#'
mvGWAS$methods(
  conduct_scan_slurm_ = function(max_num_nodes = Inf) {

    genotype_files <- list.files(path = metadata$genotype_directory, full.names = TRUE, pattern = '\\.vcf|\\.VCF')

    mean_formula = .self$metadata$mean_formula
    var_formula = .self$metadata$var_formula
    mean_null_formula = .self$metadata$mean_null_formula
    var_null_formula = .self$metadata$var_null_formula
    phenotype_df = .self$data$phenotypes

    sjob <- rslurm::slurm_apply(f = scan_vcf_file_,
                                params = dplyr::data_frame(file_name = genotype_files),
                                nodes = min(max_num_nodes, length(genotype_files)),
                                cpus_per_node = 1,
                                add_objects = list(mean_formula = mean_formula,
                                                   var_formula = var_formula,
                                                   mean_null_formula = mean_null_formula,
                                                   var_null_formula = var_null_formula,
                                                   phenotype_df = phenotype_df))


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

    usingMethods(stash_formula_, determine_keyword_use_, conduct_scan_local_, conduct_scan_slurm_)

    system <- match.arg(arg = system)

    stash_formula_(formula = mean_formula, type = 'mean')
    stash_formula_(formula = var_formula, type = 'var')

    determine_keyword_use_()

    if (!any(grepl(pattern = 'null_formula', x = names(metadata)))) {
      stop('Use MAP_gt in at least one of mean_formula or var_formula.')
    }

    switch(EXPR = system,
           'local' = conduct_scan_local_(...),
           'slurm' = conduct_scan_slurm_(...))

    return(TRUE)

  }
)



