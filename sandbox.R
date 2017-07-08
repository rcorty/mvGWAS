library(dplyr)
library(mvGWAS)

my_gwas <- mvGWAS$new(phenotype_file = 'tests/test_data/test_phenos.RDS', genotype_directory = 'tests/test_data')

my_gwas$conduct_scan(mean_formula = sbp ~ EV1 + male + DS, var_formula = ~ EV1 + male + GT, num_cores = 1)

my_gwas$qq_plot(what = 'p_LR')
my_gwas$manhattan_plot(what = 'p_LR', lowest_nlog_p_to_plot = 1)

# gwas1$apply_genomic_control(method = 'LR')
# gwas1$genomic_control_dfs
# gwas1$results %>% glimpse()
#
#
# gwas1$qq_plot(what = 'p_LR')
# gwas1$manhattan_plot(what = 'p_LR', lowest_nlog_p_to_plot = 2)



gwas1$qq_plot(what = 'p_z_DS')
gwas1$manhattan_plot(what = 'p_z_DS', lowest_nlog_p_to_plot = 1)

gwas1$apply_genomic_control(method = 'z2_DS')


