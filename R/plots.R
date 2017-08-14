#' @title Manhattan Plot
#' @name mvGWAS_manhattan_plot
#'
#' @description makes a manhattan plot for a mvGWAS
#'
#' @param what_to_plot what to plot
#' @param lowest_nlog_p_to_plot lowest neg log10 p to plot
#' @param highlight_nlog_p_above highlight neg log10 p value above this number with larger points
#' @param highlight_size how big are the highlight points
#' @param nonhighlight_alpha how transparent are the non-highlighted points
#'
#' @return the plot
#' @importFrom dplyr %>%
#'
mvGWAS$methods(

  manhattan_plot = function(what = c('p_LR', 'p_z_DS'),
                            gc,
                            lowest_nlog_p_to_plot = 3,
                            highlight_nlog_p_above = 5,
                            highlight_size = 2,
                            nonhighlight_alpha = 0.5) {

    what <- match.arg(arg = what)
    if (missing(gc)) { gc <- any(grepl(pattern = paste0(what, '.*_gc$'), x = names(results)))}

    column_regex <- paste0('^', what, '.*', if (gc) '_gc$')

    to_plot <- results %>%
      dplyr::select(dplyr::matches(column_regex), CHROM, POS) %>%
      tidyr::gather(key = test, value = p, dplyr::matches(column_regex)) %>%
      na.omit() %>%
      dplyr::mutate(nlog_p = -log(p, base = 10))
    # dplyr::mutate(test = dplyr::recode(test, mean_asymp_p = 'mean', var_asymp_p = 'var', joint_asymp_p = 'joint'))


    to_plot2 <- to_plot %>%
      dplyr::filter(nlog_p > lowest_nlog_p_to_plot, nlog_p < highlight_nlog_p_above)

    to_plot3 <- to_plot %>%
      dplyr::filter(nlog_p > highlight_nlog_p_above)

    ggplot2::ggplot(data = to_plot2, mapping = ggplot2::aes(x = POS, y = nlog_p, color = test)) +
      ggplot2::facet_grid(~CHROM, scales = 'free_x', switch = 'x', space = 'free', drop = FALSE) +
      ggplot2::geom_point(alpha = nonhighlight_alpha) +
      ggplot2::geom_point(data = to_plot3, size = highlight_size) +
      ggplot2::scale_color_manual(values = setNames(object = c('blue', 'red', 'black'),
                                                    nm = paste0(what, '_', c('mean', 'var', 'joint'), if (gc) '_gc'))) +
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
#' @description makes a qq plot for a mvGWAS
#'
#' @param what what qq to plot
#' @param log_p whether to use (negative) log10 scale or not
#' @param gc whether to use genomic controlled values or not
#'
#' @return the plot
#' @importFrom dplyr %>%
#'
mvGWAS$methods(

  qq_plot = function(what = c('p_LR', 'p_z_DS'),
                     num_quantiles = 1000,
                     top_n = 1000,
                     log_p = TRUE) {


    what <- match.arg(arg = what)

    quantiles <- seq(from = (0.5/num_quantiles), to = 1 - (0.5/num_quantiles), length.out = num_quantiles)


    d <- results %>%
      dplyr::select(dplyr::matches(what)) %>%
      tidyr::gather(key = test, value = obs_p) %>%
      dplyr::group_by(test) %>%
      # dplyr::mutate(gc = grepl(pattern = '_gc$', x = test),
      #               test = gsub(pattern = '_gc$', replacement = '', x = test),
      #               test = factor(x = test, levels = paste0(what, '_', c('mean', 'var', 'joint')))) %>%
      # dplyr::group_by(test, gc) %>%
      dplyr::do(dplyr::data_frame(exp = quantiles,
                                  obs = quantile(x = .$obs_p,
                                                 probs = quantiles,
                                                 na.rm = TRUE,
                                                 names = FALSE)))



    d %>%
      ggplot2::ggplot(mapping = ggplot2::aes(x = if (log_p) -log10(exp) else exp,
                                             y = if (log_p) -log10(obs) else obs,
                                             color = test)) +
      ggplot2::geom_abline(slope = 1, intercept = 0) +
      ggplot2::geom_point() +
      ggplot2::coord_fixed() +
      ggplot2::xlab(if (log_p) 'theoretical -log10(p)' else 'theoretical p') +
      ggplot2::ylab(if (log_p) 'observed -log10(p)' else 'observed p') +
      ggplot2::theme_minimal() +
      ggplot2::ggtitle(label = paste0(metadata$mean_alt_formula[[2]], ': QQ plot of ', what)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.ticks = ggplot2::element_blank())
  }

)
