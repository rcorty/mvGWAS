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
                                     genotype_directory = 'Yada yada yada'))

    expect_error(object = mvGWAS$new(phenotype_file = 'You yadad the best part!',
                                     genotype_directory = 'No, I mentioned the bisque.'))

  }
)


test_that(
  desc = 'Conduct a single core local mvGWAS',
  code = {

    # mean, var, and joint testing with covariates in both models
    gwas1 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas1$conduct_scan(mean_formula = sbp ~ EV1 + male + MAP_gt,
                                          var_formula = ~ EV1 + male + MAP_gt,
                                          num_cores = 1),
              class = 'mvGWAS')
    expect_is(object = gwas1$results, class = 'tbl_df')

    # mean testing with covariates in both models
    gwas2 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas2$conduct_scan(mean_formula = sbp ~ EV1 + male + MAP_gt,
                                          var_formula = ~ EV1 + male,
                                          num_cores = 1),
              class = 'mvGWAS')
    expect_is(object = gwas1$results, class = 'tbl_df')

    # var testing with covariates in both models
    gwas3 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas3$conduct_scan(mean_formula = sbp ~ EV1 + male,
                                          var_formula = ~ EV1 + male + MAP_gt,
                                          num_cores = 1),
              class = 'mvGWAS')
    expect_is(object = gwas1$results, class = 'tbl_df')

    # no testing with covariates in both models
    gwas4 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_error(object = gwas4$conduct_scan(mean_formula = sbp ~ EV1 + male,
                                             var_formula = ~ EV1 + male,
                                             num_cores = 1))
    expect_equal(object = dim(gwas4$results), c(0, 0))


    # mean, var, and joint testing with no covariates
    gwas5 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas5$conduct_scan(mean_formula = sbp ~ MAP_gt,
                                          var_formula = ~ MAP_gt,
                                          num_cores = 1),
              class = 'mvGWAS')
    expect_is(object = gwas5$results, class = 'tbl_df')

    # mean testing with no covariates
    gwas6 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas6$conduct_scan(mean_formula = sbp ~ MAP_gt,
                                          var_formula = ~ 1,
                                          num_cores = 1),
              class = 'mvGWAS')
    expect_is(object = gwas6$results, class = 'tbl_df')

    # var testing with no covariates
    gwas7 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas7$conduct_scan(mean_formula = sbp ~ 1,
                                          var_formula = ~ MAP_gt,
                                          num_cores = 1),
              class = 'mvGWAS')
    expect_is(object = gwas7$results, class = 'tbl_df')

    # no testing with no covariates
    gwas8 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_error(object = gwas8$conduct_scan(mean_formula = sbp ~ 1,
                                             var_formula = ~ 1,
                                             num_cores = 1))
    expect_equal(object = dim(gwas8$results), c(0, 0))
  }
)


test_that(
  desc = 'Conduct a multi-core local mvGWAS',
  code = {

    # mean, var, and joint testing with covariates in both models
    gwas1 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas1$conduct_scan(mean_formula = sbp ~ EV1 + male + MAP_gt,
                                          var_formula = ~ EV1 + male + MAP_gt),
              class = 'mvGWAS')
    expect_is(object = gwas1$results, class = 'tbl_df')

    # mean testing with covariates in both models
    gwas2 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas2$conduct_scan(mean_formula = sbp ~ EV1 + male + MAP_gt,
                                          var_formula = ~ EV1 + male),
              class = 'mvGWAS')
    expect_is(object = gwas1$results, class = 'tbl_df')

    # var testing with covariates in both models
    gwas3 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas3$conduct_scan(mean_formula = sbp ~ EV1 + male,
                                          var_formula = ~ EV1 + male + MAP_gt),
              class = 'mvGWAS')
    expect_is(object = gwas1$results, class = 'tbl_df')

    # no testing with covariates in both models
    gwas4 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_error(object = gwas4$conduct_scan(mean_formula = sbp ~ EV1 + male,
                                             var_formula = ~ EV1 + male))
    expect_equal(object = dim(gwas4$results), c(0, 0))


    # mean, var, and joint testing with no covariates
    gwas5 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas5$conduct_scan(mean_formula = sbp ~ MAP_gt,
                                          var_formula = ~ MAP_gt),
              class = 'mvGWAS')
    expect_is(object = gwas5$results, class = 'tbl_df')

    # mean testing with no covariates
    gwas6 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas6$conduct_scan(mean_formula = sbp ~ MAP_gt,
                                          var_formula = ~ 1),
              class = 'mvGWAS')
    expect_is(object = gwas6$results, class = 'tbl_df')


    # var testing with no covariates
    gwas7 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas7$conduct_scan(mean_formula = sbp ~ 1,
                                          var_formula = ~ MAP_gt),
              class = 'mvGWAS')
    expect_is(object = gwas7$results, class = 'tbl_df')

    # no testing with no covariates
    gwas8 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_error(object = gwas8$conduct_scan(mean_formula = sbp ~ 1,
                                             var_formula = ~ 1))
    expect_equal(object = dim(gwas8$results), c(0, 0))
  }
)


test_that(
  desc = 'Conduct an mvGWAS on a slurm cluster',
  code = {

    # mean, var, and joint testing with covariates in both models
    gwas1 <- mvGWAS$new(phenotype_file = '../test_data/test_phenos.RDS', genotype_directory = '../test_data')
    expect_is(object = gwas1$conduct_scan(mean_formula = sbp ~ EV1 + male + MAP_gt,
                                          var_formula = ~ EV1 + male + MAP_gt,
                                          system = 'slurm',
                                          max_num_nodes = 20),
              class = 'mvGWAS')
    expect_is(object = gwas1$results, class = 'tbl_df')
  }
)

