## ----example, eval=FALSE-------------------------------------------------
#  library(mvGWAS)
#  
#  my_first_gwas <- mvGWAS$new(phenotype_file = 'my_phenotype_file.csv',
#                              genotype_directory = 'genotype_dir')

## ----eval=FALSE----------------------------------------------------------
#  my_second_gwas <- mvGWAS$new(phenotype_file = 'my_phenotype_file.csv',
#                               genotype_directory = 'genotype_dir',
#                               genotype_file_pattern = 'chr17')

## ----eval=FALSE----------------------------------------------------------
#  my_gwas$conduct_scan(mean_formula = my_trait ~ PC1 + PC2 + sex + age + DS,
#                       var_formula = ~ PC1 + PC2 + sex + GT)

## ----eval=FALSE----------------------------------------------------------
#  my_gwas$conduct_scan(mean_formula = my_trait ~ PC1 + PC2 + sex + age + DS,
#                       var_formula = ~ PC1 + PC2 + sex + GT,
#                       num_cores = 2)

## ----eval=FALSE----------------------------------------------------------
#  my_gwas$qq_plot()
#  my_gwas$manhattan_plot()

## ----eval=FALSE----------------------------------------------------------
#  aaabbb
#  

