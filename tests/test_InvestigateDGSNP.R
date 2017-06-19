# library(GWASviaDGLM)
#
# file.names <- list(phen.file = '../CVD_Lange_dGWAS/1_ProjectData/2_ProcessedData/PhenotypesOfUnrelateds/JHS_data_2016-07-22_preprocessed.RDS',
#                    geno.files = paste0('../CVD_Lange_dGWAS/1_ProjectData/2_ProcessedData/DirectGenotypesOfUnrelateds/unrelateds_chr', 1:22, '.RDS'),
#                    map.file = '../CVD_Lange_dGWAS/1_ProjectData/2_ProcessedData/DirectGenotypesOfUnrelateds/genomic_map_direcly_genotyped_snps.RDS')
#
# phens <- readRDS(file = file.names[['phen.file']])
#
#
# InvestigateDGSNP(phenotype.name = 'BMI_v2',
#                  snp.name = 'rs12771692',
#                  file.names = file.names,
#                  keep.only = bpmeds_v1 == 'no')
#
#
# InvestigateDGSNP(phenotype.name = 'BMI_v3',
#                  snp.name = 'rs12621732',
#                  file.names = file.names,
#                  mean.formula = formula('BMI_v3^(-0.19) ~ EV1 + snp'),
#                  var.formula = ~ EV4 + snp,
#                  keep.only = bpmeds_v1 == 'yes')
#
# InvestigateDGSNP(phenotype.name = 'BMI_v3',
#                  snp.name = 'my.great.snp',
#                  snp.values = sample(x = 0:2, size = nrow(phens), replace = TRUE),
#                  chr = '3',
#                  file.names = file.names,
#                  modeling.df = phens)