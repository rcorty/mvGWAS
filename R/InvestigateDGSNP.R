#' #' @title InvestigateDGSNP
#' #' @name InvestigateDGSNP
#' #' @author Robert W. Corty
#' #'
#' #' @description Model and plot to investigate a directly genotyped single
#' #' nucleotide polymorphism
#' #'
#' #'
#' #' @param phenotype.name name of the phenotype
#' #' @param snp.name name of the snp
#' #' @param file.names names of phenotype and genotype files
#' #' @param mean.formula the mean formula
#' #' @param var.formula the variance formula
#' #' @param draw.ns should the n per group be plotted
#' #' @param draw.violins should the violion plot be plotted
#' #' @param draw.scatter should the scatter plot be plotted
#' #' @param draw.sd.bar should the bars indicating the standard deviation of each genotype be plotted
#' #' @param alpha.scatter alpha on the scatter plot
#' #' @param alpha.violin alpha on the violin
#' #'
#' #' @return nothing, just plots
#' #' @export
#' #'
#' InvestigateDGSNP <- function(phenotype.name,
#'                              snp.name,
#'                              snp.values = NULL,
#'                              chr = NULL,
#'                              file.names = NULL,
#'                              modeling.df = NULL,
#'                              dot.color = NULL,
#'                              keep.only = NULL,
#'                              mean.formula = NULL,
#'                              var.formula = NULL,
#'                              verbose.title = FALSE,
#'                              draw.violins = TRUE,
#'                              draw.scatter = TRUE,
#'                              draw.sd.bars = TRUE,
#'                              draw.summary.stats = TRUE,
#'                              alpha.scatter = 0.4,
#'                              alpha.violin = 0.4) {
#'
#'   if (is.null(modeling.df)) {
#'     modeling.df  <- readRDS(file = file.names[['phen.file']])
#'   }
#'
#'   if (is.null(snp.values)) {
#'     snp <- GetDGSNPByName(file.names = file.names,
#'                           snp.name = snp.name)
#'     chr <- attr(x = snp, which = 'chr')
#'     perm.for.genos <- match(modeling.df$IID, names(snp))
#'     modeling.df[['snp']] <- snp[perm.for.genos]
#'   } else {
#'     modeling.df[['snp']] <- snp.values
#'     if (is.null(chr)) {
#'       stop('If snp is provided as value (rather than name), chromosome must be provided.')
#'     }
#'   }
#'
#'
#'
#'   keep.only <- lazyeval::f_eval(f = lazyeval::f_capture(keep.only),
#'                                 data = modeling.df)
#'   if (!is.null(keep.only)) {
#'     keep.only[is.na(keep.only)] <- FALSE
#'     modeling.df <- modeling.df[keep.only, ]
#'   }
#'
#'
#'   # fit LM if mean.formula is given but var.formula is not ---------------------
#'   if (!is.null(mean.formula) & is.null(var.formula)) {
#'     stop('LM not yet implemented.')
#'   }
#'
#'   # fit DGLM if mean.formula and var.formula are both given --------------------
#'   if (!is.null(mean.formula) & !is.null(var.formula)) {
#'
#'     dglm.fit <- dglm::dglm(formula = mean.formula,
#'                            dformula = var.formula,
#'                            data = modeling.df)
#'
#'     short.resids <- residuals(object = dglm.fit)
#'     long.resids <- rep(NA, nrow(modeling.df))
#'     long.resids[as.numeric(names(short.resids))] <- short.resids
#'   }
#'
#'   plotting.df <- modeling.df
#'   if (exists('long.resids')) {
#'     plotting.df[['resid']] <- long.resids
#'   } else {
#'     plotting.df[['resid']] <- plotting.df[[as.character(phenotype.name)]]
#'   }
#'
#'   summary.stats <- plotting.df %>%
#'     dplyr::group_by(snp) %>%
#'     dplyr::summarise(mean = mean(resid, na.rm = TRUE),
#'               sd = sd(resid, na.rm = TRUE),
#'               label = paste0('n = ', sum(!is.na(resid)), '\n',
#'                              'm = ', signif(mean, 3), '\n',
#'                              's = ', signif(sd, 3)))
#'
#'   p <- ggplot2::ggplot()
#'   if (draw.violins) {
#'     p <- p + ggplot2::geom_violin(data = plotting.df,
#'                                   mapping = ggplot2::aes(x = snp,
#'                                                          y = resid,
#'                                                          group = snp),
#'                                   fill = 'lightgray',
#'                                   alpha = alpha.violin,
#'                                   color = NA)
#'   }
#'   if (draw.scatter) {
#'     p <- p +
#'       ggplot2::geom_jitter(data = plotting.df,
#'                            mapping = ggplot2::aes_string(x = 'snp', y = 'resid', col = dot.color),
#'                            width = 0.4,
#'                            height = 0,
#'                            alpha = alpha.scatter) +
#'       ggplot2::scale_alpha(guide = 'none')
#'   }
#'   if (draw.sd.bars) {
#'     p <- p +
#'       ggplot2::geom_pointrange(data = summary.stats,
#'                                mapping = ggplot2::aes(x = snp - 0.25,
#'                                              y = mean,
#'                                              ymin = mean - sd,
#'                                              ymax = mean + sd))
#'   }
#'   if (draw.summary.stats) {
#'     p <- p +
#'       ggplot2::geom_text(data = summary.stats,
#'                          mapping = ggplot2::aes(x = snp,
#'                                        y = mean,
#'                                        label = label),
#'                          hjust = 'left',
#'                          vjust = 'center',
#'                          nudge_x = 0.22)
#'   }
#'
#'
#'   if (is.null(mean.formula) & is.null(var.formula)) {
#'     title <- paste0(phenotype.name, ' by chr ', chr, ', ', snp.name)
#'   }
#'   if (!is.null(mean.formula) & is.null(var.formula)) {
#'     stop('LM not yet implemented.')
#'   }
#'   if (!is.null(mean.formula) & !is.null(var.formula)) {
#'
#'     title <- paste0('Residuals from ', Reduce(paste, deparse(mean.formula)),
#'                     '\n with residual variance ', Reduce(paste, deparse(var.formula)),
#'                     '\n with snp = chr, ', chr, ', ', snp.name)
#'
#'     if (verbose.title) {
#'
#'       all.mean.coefs <-signif(summary(dglm.fit)$coef[,c(1,4)], 3)
#'
#'       # longest.length <- max(length(as.character(all.mean.coefs)))
#'       # mean.coef.string.df <- dplyr::tbl_df(all.mean.coefs) %>%
#'       #   dplyr::mutate_all(dplyr::funs(char = stringr::str_pad(., width = longest.length))) %>%
#'       #   dplyr::select(-(1:2))
#'       #
#'       # mean.coef.string <- paste(row.names(all.mean.coefs),
#'       #                           apply(X = mean.coef.string.df,
#'       #                                 MARGIN = 1,
#'       #                                 FUN = function(r) {
#'       #                                   paste(r, collapse = '  ')
#'       #                                 }),
#'       #                           collapse = '\n')
#'       #
#'       #
#'       # title <- paste(title,
#'       #                paste(stringr::str_pad(string = c('', 'est', 'p'), width = longest.length), collapse = ''),
#'       #                mean.coef.string,
#'       #                sep = '\n')
#'
#'       title <- paste0(title,
#'                       '\n',
#'                       paste(names(all.mean.coefs[,1]), 'mean est =', all.mean.coefs[,1], collapse = ', '),
#'                       '\n',
#'                       paste(names(all.mean.coefs[,2]), 'mean p =', all.mean.coefs[,2], collapse = ', '))
#'
#'       all.var.coefs <- signif(summary(dglm.fit)$dispersion.summary$coef[,c(1,4)], 3)
#'
#'       title <- paste0(title,
#'                       '\n',
#'                       paste(names(all.var.coefs[,1]), 'mean est =', all.var.coefs[,1], collapse = ', '),
#'                       '\n',
#'                       paste(names(all.var.coefs[,2]), 'mean p =', all.var.coefs[,2], collapse = ', '))
#'
#'     }
#'     if (!verbose.title) {
#'
#'       if ('snp' %in% all.names(mean.formula)) {
#'         snp.mean.coef <- signif(summary(dglm.fit)$coef['snp',c(1,4)], 2)
#'         title <- paste0(title,
#'                         '\n', paste(c('mean est', 'p'), snp.mean.coef, sep = '=', collapse = ', '))
#'       }
#'       if ('snp' %in% all.names(var.formula)) {
#'         snp.var.coef <- signif(summary(dglm.fit)$dispersion.summary$coef['snp',c(1,4)], 2)
#'         title <- paste0(title,
#'                         '\n', paste(c('var est', 'p'), snp.var.coef, sep = '=', collapse = ', '))
#'       }
#'
#'     }
#'
#'
#'
#'
#'
#'   }
#'
#'
#'   p <- p +
#'     ggplot2::coord_cartesian(xlim = c(-0.4, 2.4)) +
#'     ggplot2::theme_bw() +
#'     ggplot2::xlab(snp.name) +
#'     ggplot2::ylab(phenotype.name) +
#'     ggplot2::ggtitle(title)
#'
#'
#'
#'
#'   return(p)
#' }