#' @title apply_genomic_control
#' @name mvGWAS_apply_genomic_control
#' @include mvGWAS.R
#'
#' @description Applies genomic control
#'
#' @return TRUE if successful, FALSE otherwise.
#'
mvGWAS$methods(
  apply_genomic_control = function(keyword = c('LR', 'DS', 'GT', 'GP'),
                                   method = c('zsquared', 'z')) {

    keyword <- match.arg(arg = keyword)
    method <- match.arg(arg = method)

    switch(EXPR = method,
           LR = apply_genomic_control_LR_(),
           z = apply_genomic_control_DS_(method = method),
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
    lambda_col_name <- paste0('gc_lambda_', test)

    lambda_df <- results %>%
      dplyr::group_by(UQ(sym(df_col_name))) %>%
      dplyr::summarise(obs_median_LR = median(x = UQ(sym(LR_col_name)), na.rm = TRUE)) %>%
      dplyr::mutate(expected_median_LR = qchisq(p = 0.5, df = UQ(sym(df_col_name))),
                    UQ(lambda_col_name) := obs_median_LR / expected_median_LR)

    LR_gc_col_name <- paste0(LR_col_name, '_gc')

    new_cols <- results %>%
      inner_join(y = lambda_df, by = df_col_name) %>%
      mutate(UQ(LR_gc_col_name) := UQ(sym(LR_col_name)) / UQ(sym(lambda_col_name))) %>%
      select(UQ(sym(LR_gc_col_name)), UQ(sym(lambda_col_name)))

    results <<- bind_cols(results, new_cols)

  }
)





#' @title apply_genomic_control_DS_
#' @name mvGWAS_apply_genomic_control
#'
#' @return TRUE if successful, FALSE otherwise.
#'
mvGWAS$methods(
  apply_genomic_control_DS_ = function(method) {

    if ('DS_z_mean' %in% names(results)) {
      gc_DS_one_test_(test = 'mean', method = method)
    }
    if ('DS_z_var' %in% names(results)) {
      gc_DS_one_test_(test = 'var', method = method)
    }

    return(TRUE)
  }
)



#' @title gc_DS_one_test_
#' @name mvGWAS_internals
#'
#' @return TRUE if successful, FALSE otherwise.
#' @importFrom dplyr %>%
#' @importFrom rlang UQ
#' @importFrom rlang sym
#'
mvGWAS$methods(
  gc_DS_one_test_ = function(test) {

    DS_z_col_name <- paste0('DS_z_', test)
    lambda_col_name <- paste0('gc_lambda_', test)
    DS_z_gc_col_name <- paste0(DS_z_col_name, '_gc')

    lambda_df <- results %>%
      dplyr::summarise(median_zsq = median(x = UQ(sym(DS_z_col_name))^2, na.rm = TRUE)) %>%
      dplyr::mutate(UQ(lambda_col_name) := median_zsq / qchisq(p = 0.5, df = 1))


    new_cols <- results %>%
      mutate(UQ(lambda_col_name) := lambda_df %>% pull(lambda_col_name),
             UQ(DS_z_gc_col_name) := UQ(sym(DS_z_col_name)) / sqrt(UQ(sym(lambda_col_name)))) %>%
      select(UQ(sym(lambda_col_name)), UQ(sym(DS_z_gc_col_name)))

    results <<- bind_cols(results, new_cols)

  }
)
