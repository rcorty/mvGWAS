#' @title manhattan_plot
#' @name mvGWAS_manhattan_plot
#' @include mvGWAS.R
#'
#' @description makes a Manhattan plot
#'
#' @return the plot
#' @importFrom dplyr %>%
#'
mvGWAS$methods(
  manhattan_plot = function(what_to_plot = c('asymp_p', 'asymp_p_gc'),
                            lowest_nlog_p_to_plot = 3,
                            highlight_nlog_p_above = 5,
                            highlight_size = 2,
                            nonhighlight_alpha = 0.5) {

    what_to_plot <- match.arg(arg = what_to_plot)

    if (what_to_plot == 'asymp_p') {

      for_plotting <- results %>%
        tidyr::gather(key = test, value = p, mean_asymp_p, var_asymp_p, joint_asymp_p) %>%
        dplyr::mutate(test = dplyr::recode(test, mean_asymp_p = 'mean', var_asymp_p = 'var', joint_asymp_p = 'joint'))

    } else {

      stop('mahnattan_plot() on gc values not yet implemented')

    }

    for_plotting2 <- for_plotting %>%
      dplyr::mutate(nlog_p = -log(p, base = 10)) %>%
      dplyr::filter(nlog_p > lowest_nlog_p_to_plot, nlog_p < highlight_nlog_p_above)

    for_plotting3 <- for_plotting %>%
      dplyr::mutate(nlog_p = -log(p, base = 10)) %>%
      dplyr::filter(nlog_p > highlight_nlog_p_above)

    ggplot2::ggplot(data = for_plotting2, mapping = ggplot2::aes(x = POS, y = nlog_p, color = test)) +
      ggplot2::facet_grid(~CHROM, scales = 'free_x', switch = 'x', space = 'free', drop = FALSE) +
      ggplot2::geom_point(alpha = nonhighlight_alpha) +
      ggplot2::geom_point(data = for_plotting3, size = highlight_size) +
      ggplot2::scale_color_manual(values = c(mean = 'blue', var = 'red', joint = 'black')) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank()) +
      ggplot2::xlab('Chromosome and position in Mb') +
      ggplot2::ylab('-log10(p)')

  }
)




#' @title qq_plot
#' @name mvGWAS_qq_plot
#'
#' @description makes a Manhattan plot
#'
#' @return the plot
#' @importFrom dplyr %>%
#'
mvGWAS$methods(
  qq_plot = function(what_to_plot,
                     theoretical = 'unif') {

    v <- results[[what_to_plot]]

    for_plotting <- dplyr::data_frame(obs = sort(v),
                                      theo = seq(from = 0, to = 1, length.out = sum(!is.na(v))))

    ggplot2::ggplot() +
      ggplot2::geom_segment(data = dplyr::data_frame(start = 0, stop = 1),
                            mapping = ggplot2::aes(x = start, y = start, xend = stop, yend = stop),
                            color = 'gray') +
      ggplot2::geom_point(data = for_plotting,
                          mapping = ggplot2::aes(x = theo, y = obs)) +
      ggplot2::ggtitle(label = paste0('QQ plot of ', what_to_plot)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                     axis.text = ggplot2::element_blank()) +
      ggplot2::xlab('Theoretical uniform distribution') +
      ggplot2::ylab('Observed p values')


  }
)