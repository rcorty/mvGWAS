context('Testing mvGWAS')

test_that(
  desc = 'Initialize an mvGWAS',
  code = {

    # proper use
    expect_is(object = mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS',
                                  genotype_directory = '../test_data'),
              class = 'mvGWAS')

    # errors
    expect_error(object = mvGWAS$new(phenotype_file = 'We went out to dinner.',
                                     genotype_directory = '../test_data'))

    expect_error(object = mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS',
                                     genotype_directory = 'Yada yada yada...'))

    expect_error(object = mvGWAS$new(phenotype_file = 'You yadad the best part!',
                                     genotype_directory = 'No, I mentioned the bisque.'))

  }
)


test_that(
  desc = 'Conduct a single core local mvGWAS',
  code = {

    # mean, var, and joint testing with covariates in both models
    gwas1 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_true(object = gwas1$conduct_scan(mean_formula = sbp ~ EV1 + male + GP,
                                            var_formula = ~ EV1 + male + GP,
                                            num_cores = 1))
    expect_is(object = gwas1$results, class = 'tbl_df')


    # mean, var, and joint testing with no covariates
    gwas2 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_true(object = gwas2$conduct_scan(mean_formula = sbp ~ EV1 + male,
                                            var_formula = ~ EV1 + male,
                                            num_cores = 1))
    expect_equal(object = dim(gwas2$results), c(0, 0))

  }
)


test_that(
  desc = 'Conduct a multi-core local mvGWAS',
  code = {

    # mean, var, and joint testing with covariates in both models
    gwas3 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_true(object = gwas3$conduct_scan(mean_formula = sbp ~ EV1 + male + GP,
                                            var_formula = ~ EV1 + male + GP))
    expect_is(object = gwas3$results, class = 'tbl_df')

    # no testing with covariates in both models
    gwas4 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_true(object = gwas4$conduct_scan(mean_formula = sbp ~ EV1 + male,
                                            var_formula = ~ EV1 + male))
    expect_equal(object = dim(gwas4$results), c(0, 0))

  }
)


test_that(
  desc = 'Conduct an mvGWAS on a slurm cluster',
  code = {

    # I have no idea why this works to test whether a slurm cluster is available, but it does
    if (as.logical(system(command = 'squeue', ignore.stdout = TRUE, ignore.stderr = TRUE))) {

      message('No slurm cluster available, so not testing submission to cluster')

    } else {

      # mean, var, and joint testing with covariates in both models
      gwas1 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
      expect_is(object = gwas1$conduct_scan(mean_formula = sbp ~ EV1 + male + DS,
                                            var_formula = ~ EV1 + male + GT,
                                            system = 'slurm'),
                class = 'mvGWAS')
      expect_is(object = gwas1$results, class = 'tbl_df')

    }
  }
)

