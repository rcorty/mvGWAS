library(GWASviaDGLM)

file.names <- list(phen.file = '../CVD_Lange_dGWAS/1_ProjectData/2_ProcessedData/PhenotypesOfUnrelateds/JHS_data_2016-07-22_preprocessed.RDS',
                   geno.files = paste0('../CVD_Lange_dGWAS/1_ProjectData/2_ProcessedData/DirectGenotypesOfUnrelateds/unrelateds_chr', 1:22, '.RDS'),
                   map.file = '../CVD_Lange_dGWAS/1_ProjectData/2_ProcessedData/DirectGenotypesOfUnrelateds/genomic_map_direcly_genotyped_snps.RDS')

phens <- readRDS(file = file.names[['phen.file']])

InvestigateDGSNP(phenotype.name = 'sbp_v1_adj',
                 snp.name = 'rs11143868',
                 file.names = file.names,
                 mean.formula = sbp_v1_adj ~ EV1 + snp,
                 var.formula = ~ EV4 + snp,
                 keep.only = sex == 'male',
                 verbose.title = TRUE)