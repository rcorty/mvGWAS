library(dplyr)
library(mvGWAS)

gwas1 <- mvGWAS$new(phenotype_file = 'tests/test_data/test_phenos.RDS', genotype_directory = 'tests/test_data')

gwas1$conduct_scan(mean_formula = sbp ~ EV1 + male + DS, var_formula = ~ EV1 + male + DS)

gwas1$qq_plot(what = 'p_LR')
gwas1$manhattan_plot(what = 'p_LR', lowest_nlog_p_to_plot = 2)

gwas1$apply_genomic_control(keyword = 'LR')
gwas1$genomic_control_dfs
gwas1$results %>% glimpse()


gwas1$qq_plot(what = 'p_LR')
gwas1$manhattan_plot(what = 'p_LR', lowest_nlog_p_to_plot = 2)



gwas1$qq_plot(what = 'p_z_DS')
gwas1$manhattan_plot(what = 'p_z_DS', lowest_nlog_p_to_plot = 1)

gwas1$apply_genomic_control(keyword = 'z_DS')



# gwas1$apply_genomic_control(keyword = 'DS')
# gwas1$genomic_control_dfs
# gwas1$results %>% glimpse()
#
#
