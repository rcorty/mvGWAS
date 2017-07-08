#' @title apply_genomic_control
#' @name mvGWAS_apply_genomic_control
#' @include mvGWAS.R
#'
#' @description Applies genomic control
#'
#' @return TRUE if successful, FALSE otherwise.
#'
mvGWAS$methods(
  apply_genomic_control = function(method = c('LR', 'z_DS', 'z2_DS')) {

    method <- match.arg(arg = method)

    switch(EXPR  = method,
           LR    = apply_genomic_control_LR_(),
           z_DS  = apply_genomic_control_z_DS_(),
           z2_DS = apply_genomic_control_z2_DS_(),
           FALSE)
  }
)






#' @title apply_genomic_control_LR_
#' @name mvGWAS_apply_genomic_control
#'
#' @return TRUE if successful, FALSE otherwise.
#'
mvGWAS$methods(
  apply_genomic_control_LR_ = function() {

    if ('LR_mean' %in% names(results)) {
      gc_LR_one_test_(test = 'mean')
    }
    if ('LR_var' %in% names(results)) {
      gc_LR_one_test_(test = 'var')
    }
    if ('LR_joint' %in% names(results)) {
      gc_LR_one_test_(test = 'joint')
    }

    return(TRUE)
  }
)




#' @title gc_LR_one_test_
#' @name mvGWAS_internals
#'
#' @return TRUE if successful, FALSE otherwise.
#' @importFrom dplyr %>%
#' @importFrom rlang UQ
#' @importFrom rlang sym
#'
mvGWAS$methods(
  gc_LR_one_test_ = function(test) {

    LR_col_name <- paste0('LR_', test)
    df_col_name <- paste0('df_', test)

    LR_col_name_gc <- paste0(LR_col_name, '_gc')
    p_LR_col_name_gc <- paste0('p_', LR_col_name_gc)


    lambda_df <- results %>%
      dplyr::group_by(UQ(sym(df_col_name))) %>%
      dplyr::summarise(obs_median_LR = median(x = UQ(sym(LR_col_name)), na.rm = TRUE)) %>%
      dplyr::mutate(expected_median_LR = qchisq(p = 0.5, df = UQ(sym(df_col_name))),
                    lambda = obs_median_LR / expected_median_LR)

    results <<- results %>%
      inner_join(y = lambda_df, by = df_col_name) %>%
      mutate(UQ(LR_col_name_gc) := UQ(sym(LR_col_name)) / lambda,
             UQ(p_LR_col_name_gc) := pchisq(q = LR_gc,  df = UQ(sym(df_col_name)),  lower.tail = FALSE))

    genomic_control_dfs[[paste0('LR_', test)]] <<- na.omit(lambda_df)

  }
)


#' #' @title apply_genomic_control_DS_
#' #' @name mvGWAS_apply_genomic_control
#' #'
#' #' @return TRUE if successful, FALSE otherwise.
#' #'
# mvGWAS$methods(
#   apply_genomic_control_z_DS_ = function(method) {
#
#     if ('z_DS_mean' %in% names(results)) {
#       gc_z_DS_one_test_(test = 'mean', method = method)
#     }
#     if ('z_DS_var' %in% names(results)) {
#       gc_z_DS_one_test_(test = 'var', method = method)
#     }
#
#     return(TRUE)
#   }
# )
#
#
#
#' #' @title gc_DS_one_test_
#' #' @name mvGWAS_internals
#' #'
#' #' @return TRUE if successful, FALSE otherwise.
#' #' @importFrom dplyr %>%
#' #' @importFrom rlang UQ
#' #' @importFrom rlang sym
#' #'
# mvGWAS$methods(
#   gc_z_DS_one_test_ = function(test) {
#
#     z_DS_col_name <- paste0('z_DS_', test)
#
#     z_DS_col_name_gc <- paste0(z_DS_col_name, '_gcz')
#     p_z_DS_col_name_gc <- paste0('p_', z_DS_z_col_name_gc)
#
#     quantiles <- results %>% dplyr::pull(UQ(sym(z_DS_col_name))^2) %>% quantile(probs = c(0.4, 0.5, 0.6))
#
#     results <<- results %>%
#       mutate(UQ(z_DS_col_name_gcz))
#       mutate(z_DS_gc := UQ(sym(z_DS_col_name)) / sqrt(lambda_df$lambda[1])) %>%
#       pull(z_DS_gc)
#
#     genomic_control_dfs[[paste0('z_DS_', test)]] <<- dplyr::data_frame(lambda = lambda)
#
#   }
# )


#' @title apply_genomic_control_z2_DS_
#' @name mvGWAS_apply_genomic_control
#'
#' @return TRUE if successful, FALSE otherwise.
#'
mvGWAS$methods(
  apply_genomic_control_z2_DS_ = function() {

    if ('z_DS_mean' %in% names(results)) {
      gc_z2_DS_one_test_(test = 'mean')
    }
    if ('z_DS_var' %in% names(results)) {
      gc_z2_DS_one_test_(test = 'var')
    }

    return(TRUE)
  }
)




#' @title gc_LR_one_test_
#' @name mvGWAS_internals
#'
#' @return TRUE if successful, FALSE otherwise.
#' @importFrom dplyr %>%
#' @importFrom rlang UQ
#' @importFrom rlang sym
#'
mvGWAS$methods(
  gc_z2_DS_one_test_ = function(test) {

    z_DS_col_name <- paste0('z_DS_', test)

    z_DS_col_name_gc <- paste0(z_DS_col_name, '_gcz2')
    p_z_DS_col_name_gc <- paste0('p_', z_DS_col_name_gc)


    lambda <- results %>%
      pull(UQ(sym(z_DS_col_name))) %>%
      `^`(2) %>%
      median(na.rm = TRUE) %>%
      `/`(qchisq(p = 0.5, df = 1))

    results <<- results %>%
      mutate(UQ(z_DS_col_name_gc) := UQ(sym(z_DS_col_name)) / lambda,
             UQ(p_z_DS_col_name_gc) := dplyr::if_else(condition = UQ(sym(z_DS_col_name_gc)) < 0,
                                                     true = pnorm(q = UQ(sym(z_DS_col_name_gc))),
                                                     false = pnorm(q = UQ(sym(z_DS_col_name_gc)), lower.tail = TRUE),
                                                     missing = NA_real_))


    genomic_control_dfs[[paste0('z2_DS_', test)]] <<- dplyr::data_frame(lambda = lambda)

  }
)
