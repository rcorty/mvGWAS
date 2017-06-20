#' #' @title QQPlot Function
#' #' @name QQPlot
#' #' @author Robert Corty
#' #'
#' #' @description This function takes a vector of values, sorts them, and plots them
#' #' against the expected distribution specified by the \code{theoretical} argument.
#' #' As \code{v} gets longer, the less interesting end(s) of the distribution get
#' #' sparisifed to save plotting time.
#' #'
#' #' @param v values to be plotted
#' #' @param theoretical theoretical distribution of the values \code{v},
#' #' defaults to 'unif'.  Other possible values are 'normal'
#' #' @param plot.neg.log.ten Whether \code{v} and its theoretical values should be
#' #' negative-log-tenned before plotting.
#' #'
#' #' @details For uniform theoretical distribution, the 'less interesting' end
#' #' is considered the end near 1.  The values are always sorted at the start.
#' #' The 'most interesting' 1e3 values are always plotted.  The values up through
#' #' 1e4 are plotted in a '1 out of every 10' fashion.  The values up through 1e5
#' #' are plotted in a '1 out of every 100 fashion'.  The values up through 1e6 are
#' #' plotted in a '1 out of every 1000 fashion'.  The values up through 1e7 are
#' #' plotted in a '1 out of every 10000' fashion.
#' #'
#' #' @return No return.  Just the plot.
#' #'
#' #' @export
#' #'
#' QQPlot <- function(v,
#'                    theoretical = 'unif',
#'                    do.neg.log.ten = TRUE,
#'                    color = 'black') {
#'
#'   n <- length(v) - sum(is.na(v))
#'
#'   if (length(v) - sum(is.na(v)) > 1e7) {
#'     stop('QQPlot can only handle vectors of size less than 1e7 at this point.')
#'   }
#'
#'   if (n < 1e3) {
#'     show <- 1:n
#'   } else {
#'     if (n < 1e4) {
#'       show <- c(seq(from = 1, to = 1e3, by = 1e0),
#'                 seq(from = 1e3 + 1, to = n, by = 1e1))
#'     } else {
#'       if (n < 1e5) {
#'         show <- c(seq(from = 1, to = 1e3, by = 1e0),
#'                   seq(from = 1e3 + 1, to = 1e4, by = 1e1),
#'                   seq(from = 1e4 + 1, to = n, by = 1e2))
#'       } else {
#'         if (n < 1e6) {
#'           show <- c(seq(from = 1, to = 1e3, by = 1e0),
#'                     seq(from = 1e3 + 1, to = 1e4, by = 1e1),
#'                     seq(from = 1e4 + 1, to = 1e5, by = 1e2),
#'                     seq(from = 1e5 + 1, to = n, by = 1e3))
#'         } else {
#'           if (n < 1e7) {
#'             show <- c(seq(from = 1, to = 1e3, by = 1e0),
#'                       seq(from = 1e3 + 1, to = 1e4, by = 1e1),
#'                       seq(from = 1e4 + 1, to = 1e5, by = 1e2),
#'                       seq(from = 1e5 + 1, to = 1e6, by = 1e3),
#'                       seq(from = 1e6 + 1, to = n, by = 1e4))
#'           } else {
#'             stop('QQPlot can only take 1e7 points')
#'           }
#'         }
#'       }
#'     }
#'   }
#'
#'   if (theoretical == 'unif') {
#'
#'     theo <- -log10(seq(from = 0, to = 1, length.out = n))
#'
#'     if (do.neg.log.ten) {
#'       v <- sort(x = -log10(v), na.last = NA, decreasing = TRUE)
#'     } else {
#'       v <- sort(x = v, na.last = NA, decreasing = TRUE)
#'     }
#'
#'     data.frame(emp = v[show], theo = theo[show]) %>%
#'       ggplot2::ggplot(ggplot2::aes(x = theo, y = emp)) +
#'       ggplot2::geom_point(color = color) +
#'       ggplot2::geom_abline(slope = 1, intercept = 0) +
#'       ggplot2::coord_fixed(ratio = 1) +
#'       ggplot2::theme_bw()
#'
#'   } else {
#'     stop('Theoretical distributions other than uniform not yet implemented.')
#'   }
#' }
