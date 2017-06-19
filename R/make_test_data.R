# library(tidyverse)
# library(vcfR)
#
#
#
# x <- read.vcfR(file = '../genotypes/chr22_chunk_3_mac_gt_30.vcf')
#
#
# write.vcf(x = x[0001:0200, 1:101], file = 'test_genos1.vcf')
# write.vcf(x = x[1001:1200, 1:101], file = 'test_genos2.vcf')
# write.vcf(x = x[2001:2200, 1:101], file = 'test_genos3.vcf')
# write.vcf(x = x[3001:3200, 1:101], file = 'test_genos4.vcf')
# write.vcf(x = x[4001:4200, 1:101], file = 'test_genos5.vcf')
# write.vcf(x = x[5001:5200, 1:101], file = 'test_genos6.vcf')
# write.vcf(x = x[6001:6200, 1:101], file = 'test_genos7.vcf')
#
#
# ids <- colnames(x@gt)[2:101]
#
# d <- data_frame(ID = ids,
#                 sbp = rnorm(n = 100),
#                 EV1 = rnorm(n = 100),
#                 male = sample(x = c(0, 1), size = 100, replace = TRUE))
#
# saveRDS(object = d,
#         file = 'test_phenos.RDS')
