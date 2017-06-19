#' @title ManhattanPlot
#' @name ManhattanPlot
#' @author Robert Corty
#'
#' @param df The tbl_df to be plotted.  Must have columns called 'chr', 'pos',
#' 'mean.p', and 'var.p'.
#' @param lower.bound The lowest value to be plotted.  This number will be
#' interpreted on the negative log ten scale, so 2 implies that the least
#' significant p-value to be plotted will be 0.01.
#' @param hit.cutoff The lowest value that gets a bigger dot, drawing
#' attention to it.  Default to \code{-log10(0.05/nrow(df))}.
#' @param hit.alpha Alpha value (opacity) for SNPs above the hit cutoff.
#' Defaults to 0.8
#' @param non.hit.alpha Alpha value (opacity) for SNPs below the hit cutoff.
#' Defaults to 0.3.
#' @param hit.point.size Size of points above the hit cutoff.  Defaults
#' to 1.
#' @param non.hit.point.size Size of poitns below the hit cutoff.  Defaults
#' to 0.3.
#' @param do.neg.log.ten Whether this ploting function should compute the
#' -log10 of the p-values supplied in \code{df}.  Defaults to TRUE, which
#' is valid when p-values are supplied.  -log10 of p-values may be provided
#' in some cases, e.g. when using this plotting function in the context of
#' a Shiny app that uses brushing to select a SNP.
#' @param use.plotly should plot be wrapped in a plotly object?
#'
#' @return a plot
#' @export
#'
ManhattanPlot <- function(df,
                          lower.bound = 2,
                          hit.cutoff = -log10(0.05/nrow(df)),
                          hit.alpha = 0.8,
                          non.hit.alpha = 0.3,
                          hit.point.size = 1,
                          non.hit.point.size = 0.5,
                          do.neg.log.ten = TRUE,
                          use.plotly = FALSE) {

  # Make sure all the right columns are present
  # May want to also check their types
  stopifnot(all(c('chr', 'pos', 'mean.p', 'var.p') %in% names(df)))

  # Make sure no p-values are infinite
  stopifnot(!any(is.infinite(df$mean.p), is.infinite(df$var.p)))

  # switch to neg log ten scale if needed
  if (do.neg.log.ten) {
    df$mean.p <- -log10(df$mean.p)
    df$var.p <- -log10(df$var.p)
  }

  p <- ggplot(data = df, mapping = aes(x = pos)) +
    facet_grid(~chr, scales = 'free_x', switch = 'x', space = 'free', drop = FALSE) +
    geom_point(data = subset(df, mean.p > lower.bound & !is.na(mean.p)), mapping = aes(y = mean.p), color = 'blue', size = non.hit.point.size, alpha = non.hit.alpha) +
    geom_point(data = subset(df, mean.p > hit.cutoff & !is.na(mean.p)), mapping = aes(y = mean.p), color = 'blue', size = hit.point.size, shape = 21, alpha = hit.alpha) +
    geom_point(data = subset(df, var.p > lower.bound & !is.na(var.p)), mapping = aes(y = var.p), color = 'red', size = non.hit.point.size, alpha = non.hit.alpha) +
    geom_point(data = subset(df, var.p > hit.cutoff & !is.na(var.p)), mapping = aes(y = var.p), color = 'red', size = hit.point.size, shape = 21, alpha = hit.alpha) +
    scale_y_continuous(expand = c(0.02, 0)) +
    ggthemes::theme_hc() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    xlab('Chromosome and position in Mb') +
    ylab('-log10(p)')

  if (use.plotly) {
    return(plotly::ggplotly(p))
  } else {
    return(p)
  }
}