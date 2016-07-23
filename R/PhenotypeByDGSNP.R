#' @title PhenotypeByDGSNP Function
#' @name PhenotypeByDGSNP
#' @author Robert Corty
#'
#' @description This function plots the specified phenotype
#' stratified by the genotype at the specified SNP.
#' options [a, b, and c] must have been set a priori.
#'
#' @param phenotype.name the phenotype to plot
#' @param snp.name the rsid of the SNP to stratify by
#' @param draw.ns should the plot show the number of individuals in each
#' genotype group?  Defaults to FALSE
#'
#' @return No return.  Just the plot.
#'
#' @export
#'
PhenotypeByDGSNP <- function(phenotype.name,
                             snp.name,
                             file.names,
                             box.cox.lambda = 1,
                             mean.formula = NULL,
                             var.formula = NULL,
                             draw.ns = FALSE,
                             violin = FALSE,
                             scatter = TRUE,
                             sd.bars = TRUE,
                             alpha.scatter = 0.4,
                             alpha.violin = 0.4,
                             ...) {

  # todo: figure out how to setup interface to deal with specification of
  # response via phenotype.name and box.cox.lambda or via
  # second element of mean.formula, e.g. mean.formula[[2]]

  snp <- GetDGSNPByName(file.names = file.names,
                        snp.name = snp.name)
  chr <- attr(x = snp, which = 'chr')

  phenos <- readRDS(file = file.names[['phen.file']])

  perm.for.genos <- match(phenos$IID, names(snp))

  snp <- snp[perm.for.genos]

  df <- data_frame(snp = snp,
                   phen = if(box.cox.lambda == 0) {
                     log(phenos[[phenotype.name]])
                   } else {
                     phenos[[phenotype.name]]^box.cox.lambda
                   })

  dglm.fit <- NULL
  if (all(!is.null(mean.formula), !is.null(var.formula))) {
    dglm.fit <- dglm(formula = mean.formula,
                     dformula = var.formula,
                     data = phenos %>% mutate(snp = snp),
                     ...)
  }


  # ink.calculator <- if (ink.prop.to.log.sample.size) {
  #   ifelse(count == 1,
  #          yes = 1,
  #          no = 1/(1 + log10(count)))
  # } else {
  #   0.4
  # }

  recap <- df %>%
    group_by(snp) %>%
    summarise(count = n(),
              phen.mean = mean(phen, na.rm = TRUE),
              phen.sd = sd(phen, na.rm = TRUE),
              message = paste0('n = ', count, '\n',
                               'm = ', signif(phen.mean, 3), '\n',
                               's = ', signif(phen.sd, 3)))


  df <- recap %>% left_join(df, by = 'snp')

  p <- ggplot()
  if (violin) {
    p <- p + geom_violin(data = df,
                         mapping = aes(x = snp, y = phen, group = snp),
                         fill = 'lightgray',
                         alpha = alpha.violin,
                         color = NA)
  }
  if (scatter) {
    p <- p +
      geom_jitter(data = df, mapping = aes(x = snp, y = phen),
                  width = 0.4, height = 0, alpha = alpha.scatter) +
      scale_alpha(guide = 'none')
  }
  if (sd.bars) {
    p <- p +
      geom_pointrange(data = recap,
                      mapping = aes(x = snp - 0.25,
                                    y = phen.mean,
                                    ymin = phen.mean - phen.sd,
                                    ymax = phen.mean + phen.sd))
  }
  if (draw.ns) {
    p <- p +
      geom_text(data = recap,
                mapping = aes(x = snp,
                              y = phen.mean,
                              label = message),
                hjust = 'left',
                vjust = 'center',
                nudge_x = 0.22)
  }

  if (is.null(mean.formula)) {
    title <- paste0(phenotype.name, '^', box.cox.lambda, ' by chr ', chr, ', ', snp.name)
  } else {
    mean.snp.coef <- signif(summary(dglm.fit)$coef['snp',c(1,2,4)], 2)
    var.snp.coef <- signif(summary(dglm.fit)$dispersion.summary$coef['snp',c(1,2,4)], 2)
    title <- paste0(phenotype.name, '^', box.cox.lambda, ' by chr ', chr, ', ', snp.name,
                    '\n', paste(c('mean est', 'se', 'p'), mean.snp.coef, sep = '=', collapse = ', '),
                    '\n', paste(c('var est', 'se', 'p'), var.snp.coef, sep = '=', collapse = ', '))
  }

  p <- p +
    coord_cartesian(xlim = c(-0.3, 2.3)) +
    theme_bw() +
    xlab(snp.name) +
    ylab(phenotype.name) +
    ggtitle(title)

  return(p)
}