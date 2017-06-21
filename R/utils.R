#' @title tryNA
#' @description silently get what you can out of an expression
#'
#' @param expr the expression to try to evaluate
#'
#' @return If the expression evaluates without error, the evaluated expression.
#' If there is an error in evaluating the expression, NA.
#' Warnings are suppresssed.
#'
#' @export
#'
tryNA <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NA,
                            finally = NA))
}


#' @title tryNULL
#' @description silently get what you can out of an expression
#'
#' @param expr the expression to try to evaluate
#'
#' @return If the expression evaluates without error, the evaluated expression.
#' If there is an error in evaluating the expression, NULL.
#' Warnings are suppresssed.
#'
#' @export
#'
tryNULL <- function(expr) {
  suppressWarnings(tryCatch(expr = expr,
                            error = function(e) NULL,
                            finally = NULL))
}



clarify_gts <- function(gts, drop_gts_w_fewer_than_x_obs) {

  # if hets in both 'directions' are present, collapse them to one
  # may need to deal with '0\1' or '0/1' at some point
  if (all('0|1' %in% gts, '1|0' %in% gts_raw)) {
    gts <- replace(x = gts, list = gts == '1|0', values = '0|1')
  }

  # if there's any very rare GT, drop it
  if (any(table(gts) < drop_gts_w_fewer_than_x_obs)) {
    bad_gts <- names(table(gts))[table(gts) < drop_gts_w_fewer_than_x_obs]
    gts <- replace(x = gts, list = which(gts == bad_gts), values = NA)
  }

  # if there's only one level, move on
  # have to check again after possibly dropping a rare GT
  if (length(unique(gts)) == 1) { return(NULL) }

  # after all that checking and pruning, finally turn gts into factor
  gts <- factor(gts)
}