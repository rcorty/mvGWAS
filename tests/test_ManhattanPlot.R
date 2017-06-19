# library(dplyr)
# library(ggplot2)
# library(GWASviaDGLM)
#
# n.chr <- 5
# n.snps.per.chr <- 50
# n.snps <- n.chr * n.snps.per.chr
#
# fake.result <- data_frame(chr = rep(c(1:n.chr), each = n.snps.per.chr),
#                           pos = rep(1:n.snps.per.chr, times = n.chr) + rnorm(n = n.snps, sd = 0.3),
#                           mean.p = runif(n = n.snps),
#                           var.p = runif(n = n.snps),
#                           snp.name = sample(x = LETTERS, size = n.snps, replace = TRUE))
#
#
# ManhattanPlot(df = fake.result, lower.bound = 0, hit.cutoff = 1.5)
#
# ManhattanPlot(df = fake.result %>% mutate(mean.p = -log10(mean.p),
#                                           var.p = -log10(var.p)),
#               do.neg.log.ten = FALSE,
#               lower.bound = 0,
#               hit.cutoff = 1.5)
#
# ManhattanPlot(df = fake.result, lower.bound = 0, hit.cutoff = 1.5, use.plotly = TRUE)
