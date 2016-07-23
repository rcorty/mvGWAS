#' @title ManhattanPlot
#' @name ManhattanPlot
#' @author Robert Corty
#'
#' @param df
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
ManhattanPlot <- function(df,
                          lower.bound = 2,
                          hit.cutoff = 4,
                          hit.alpha = 0.8,
                          non.hit.alpha = 0.3,
                          hit.point.size = 1,
                          non.hit.point.size = 0.5,
                          ps.already.logged = FALSE) {

  if (!all(c('chr', 'pos', 'mean.p', 'var.p') %in% names(df))) {
    stop("tbl_df to make Manhattan Plot must contain columns 'chr', 'pos', 'mean.p', and 'var.p'")
  }

  df$mean.p <- -log10(df$mean.p)
  df$var.p <- -log10(df$var.p)

  ggplot(data = df, mapping = aes(x = pos)) +
    facet_grid(~chr, scales = 'free_x', switch = 'x', space = 'free') +
    geom_point(data = subset(df, mean.p > lower.bound & !is.na(mean.p)), mapping = aes(y = mean.p), color = 'blue', size = non.hit.point.size, alpha = non.hit.alpha) +
    geom_point(data = subset(df, mean.p > hit.cutoff & !is.na(mean.p)), mapping = aes(y = mean.p), color = 'blue', size = hit.point.size, shape = 21, alpha = hit.alpha) +
    geom_point(data = subset(df, var.p > lower.bound & !is.na(var.p)), mapping = aes(y = var.p), color = 'red', size = non.hit.point.size, alpha = non.hit.alpha) +
    geom_point(data = subset(df, var.p > hit.cutoff & !is.na(var.p)), mapping = aes(y = var.p), color = 'red', size = hit.point.size, shape = 21, alpha = hit.alpha) +
    scale_y_continuous(expand = c(0.01, 0)) +
    ggthemes::theme_hc() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()) +
    xlab('Chromosome and position in Mb') +
    ylab('-log10(p)')
}