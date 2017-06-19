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
