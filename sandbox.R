library(mvGWAS)

gwas1 <- mvGWAS$new(phenotype_file = 'tests/test_data/test_phenos.RDS', genotype_directory = 'tests/test_data')

gwas1$conduct_scan(mean_formula = sbp ~ EV1 + male + DS,
                   var_formula = ~ EV1 + male + GT)

gwas1$results %>% glimpse()

gwas1$apply_genomic_control(keyword = 'DS')



