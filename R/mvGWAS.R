#' @title Reference Class for Mean-Variance GWASs
#' @name mvGWAS
#'
#' @description A reference class implementing a Mean-Variance Genome Wide Association Study
#'
#' @field created_at when the mvGWAS was created
#' @field metadata information about the mvGWAS
#' @field data the data used to conduct the mvGWAS
#' @field null_model the model fit without any SNPs
#' @field results the results of the mvGWAS
#' @field genomic_control_dfs data_frames holding the details of genomic control corrections applied
#'
#' @export mvGWAS
#' @exportClass mvGWAS
#'
mvGWAS <- setRefClass(
  Class = 'mvGWAS',
  fields = list(created_at = 'character',
                metadata = 'list',
                data = 'list',
                null_model = 'ANY',    # todo: make this DGLM, but need to import it beforehand or something
                results = 'data.frame',
                genomic_control_dfs = 'list'))


# could poentially lock more fields, but need to separate out things that can change from those that can't
mvGWAS$lock('created_at')


#' @title initialize
#' @name mvGWAS
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
                        genotype_directory,
                        genotype_file_pattern = '.*\\.(vcf|VCF)') {

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
    # definitely a dir containing at least one .VCF or .VCF.GZ file
    # maybe it can contain only that type of files?

    created_at <<- as.character(Sys.time())

    metadata <<- list(phenotype_file = phenotype_file,
                      genotype_directory = genotype_directory,
                      genotype_file_pattern = genotype_file_pattern,
                      available_keywords = c('GT', 'DS', 'GP'))

    data <<- list(phenotypes = phenotypes)

  }
)



# hack to make R CMD CHECK happy
mvGWAS$methods(
  new = initialize
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






