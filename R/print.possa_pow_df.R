#'@title Print pow results
#'
#'@description Prints, in a readable manner, the main information from the list
#'  created by the \code{\link[POSSA:pow]{POSSA::pow}} function, calling
#'  \code{\link{print.possa_pow_df}} for each of the POSSA power information
#'  data frames in the list. This is an extension (method) of the base R
#'  \code{\link{print}} function, so it can be called simply as \code{print()}.
#'@param x The \code{\link{list}} returned by the
#'  \code{\link[POSSA:pow]{POSSA::pow}} function.
#'@param round_to Number of fractional digits to round to, for the displayed
#'  numbers. The default is the value passed from the
#'  \code{\link[POSSA:pow]{POSSA::pow}} function (stored in the returned list).
#'@param ... (Allow additional arguments for technical reasons.)
#'
#'@return Returns nothing (\code{NULL}); this method serves only to print
#'  information to the console.
#'
#'@seealso \code{\link{pow}}, \code{\link{print.possa_pow_df}}
#'
#' @export
print.possa_pow_list = function(x,
                                round_to = NA,
                                ...) {
    if (is.na(round_to)) {
        round_to = x$arguments$round_to
    }
    pow_df_names = c()
    for (el_name in names(x)) {
        if ('possa_pow_df' %in% class(x[[el_name]])) {
            pow_df_names = c(pow_df_names, el_name)
        }
    }
    cat('\033[0;34m# POSSA pow() results #\033[0m', fill = TRUE)
    for (possa_pow_el_name in pow_df_names) {
        possa_pow_list_elem = x[[possa_pow_el_name]]
        if (length(pow_df_names) > 1) {
            cat(
                'GROUP: \033[0;40m',
                possa_pow_el_name,
                '\033[0m',
                fill = TRUE,
                sep = ''
            )
        }
        print.possa_pow_df(possa_pow_list_elem,
                           round_to = round_to,
                           possa_title = FALSE,
                           ...)
    }
}

#'@title Print pow results data frame
#'
#'@description Prints, in a readable manner, the main information from any of
#'  the data frames containing power information from the list created by the
#'  \code{\link[POSSA:pow]{POSSA::pow}} function. This is an extension (method)
#'  of the base R \code{\link{print}} function, so it can be called simply as
#'  \code{print()}.
#'@param x Power information \code{\link{data.frame}} included in the list
#'  returned by the \code{\link[POSSA:pow]{POSSA::pow}} function.
#'@param round_to Number of fractional digits to round to, for the displayed
#'  numbers (\code{5} by default).
#'@param possa_title Set to \code{FALSE} to omit title printing.
#'@param ... (Allow additional arguments for technical reasons.)
#'
#'@return Returns nothing (\code{NULL}); this method serves only to print
#'  information to the console.
#'
#'@seealso \code{\link{pow}}, \code{\link{print.possa_pow_list}}
#'
#' @export
print.possa_pow_df = function(x,
                              round_to = 5,
                              possa_title = TRUE,
                              ...) {
    cat(format.possa_pow_df(
        x = x,
        round_to = round_to,
        possa_title = possa_title,
        ...
    ))
}

format.possa_pow_df = function(x,
                               round_to = 5,
                               possa_title = TRUE,
                               ...) {
    if (possa_title) {
        cat('\033[0;34m# POSSA pow() results (df) #\033[0m', fill = TRUE)
    }
    df_pow = x
    df_nrow = nrow(df_pow)
    mlook = df_nrow - 1
    looks = 1:mlook
    if (!(all(names(df_pow)[1:3] == c('look', 'n', 'n_rate')) &&
          df_pow$look[df_nrow] == 'totals' &&
          all(df_pow$look[1:(df_nrow - 1)] == 1:(df_nrow - 1)))) {
        warning(
            'Caution!: This data.frame does not seem to have the proper structure',
            ' as created via POSSA::sim().'
        )
    }
    df_p_names = c()
    df_p_nams_all = c()
    for (colnam in names(df_pow)[startsWith(names(df_pow), 'alpha_local_')]) {
        if ('possa_p' %in% class(df_pow[[colnam]])) {
            df_p_names = c(df_p_names, colnam)
            df_p_nams_all = c(df_p_nams_all, colnam)
        } else if ('possa_p_nonstopper' %in% class(df_pow[[colnam]])) {
            df_p_nams_all = c(df_p_nams_all, colnam)
        }
    }
    df_p_names = gsub('alpha_local_', '', df_p_names)
    df_p_nams_all = gsub('alpha_local_', '', df_p_nams_all)
    df_fut_names = c()
    for (futcol in names(df_pow)[startsWith(names(df_pow), 'futil_local_')]) {
        if ('possa_futility' %in%
            class(df_pow[[futcol]])) {
            df_fut_names = c(df_fut_names, futcol)
        }
    }
    cat(
        'N(average-total) = \033[0;4m',
        ro(df_pow$n_avg_prop_0[df_nrow], 1, leading_zero = TRUE),
        '\033[0m (if H0 true) or \033[0;4m',
        ro(df_pow$n_avg_prop_1[df_nrow], 1, leading_zero = TRUE),
        '\033[0m (if H1 true)',
        sep = '',
        fill = TRUE
    )
    fut_text = ''
    for (p_nam in df_p_nams_all) {
        if (p_nam %in% df_p_names) {
            nonstp = ''
            l_a_descr = '\nLocal alphas: '
        } else {
            nonstp = '\033[0;40;3mnon-stopper: \033[0;40;0m'
            l_a_descr = '\nLocal alphas (secondary): '
        }
        p_a_locals = ro(df_pow[[paste0('alpha_local_', p_nam)]][looks], round_to)
        p_a_locals = ifelse(p_a_locals == '0',
                            'none',
                            paste0('\033[0;1m', p_a_locals, '\033[0m'))
        if (mlook > 1) {
            sign_ratios = paste0(
                '\n\033[0;3mLikelihoods of significance if H0 true: ',
                paste(paste0(
                    '(',
                    looks,
                    ') ',
                    ro(df_pow[[paste0('ratio_sign_', p_nam, '_h0')]][looks], round_to),
                    collapse = '; '
                )),
                '\nLikelihoods of significance if H1 true: ',
                paste(paste0(
                    '(',
                    looks,
                    ') ',
                    ro(df_pow[[paste0('ratio_sign_', p_nam, '_h1')]][looks], round_to),
                    collapse = '; '
                )),
                '\033[0m'
            )
        } else {
            sign_ratios = ''
        }
        toprint = paste0(
            '(',
            nonstp,
            '\033[0;40m',
            p_nam,
            '\033[0m) Type I error: \033[1;31m',
            ro(df_pow[[paste0('ratio_sign_', p_nam, '_h0')]][df_nrow], round_to),
            '\033[0m; Power: \033[0;32m',
            ro(df_pow[[paste0('ratio_sign_', p_nam, '_h1')]][df_nrow], round_to),
            '\033[0m',
            l_a_descr,
            paste(paste0(
                '(',
                looks,
                ') ',
                p_a_locals,
                collapse = '; '
            )),
            sign_ratios
        )
        cat(toprint, fill = TRUE)
        if (paste0('futil_local_',
                   p_nam) %in% df_fut_names) {
            df_pow[[paste0('futil_local_',
                           p_nam)]][is.na(df_pow[[paste0('futil_local_', p_nam)]])] = 1
            df_pow[[futcol]][is.na(df_pow[[futcol]])] = 1
            if (all(df_pow[[paste0('futil_local_', p_nam)]][looks[-mlook]] == 1)) {
                cat('Futility bounds: none', fill = TRUE)
            } else {
                p_fa_locals = ro(df_pow[[paste0('futil_local_', p_nam)]][looks[-mlook]], round_to)
                p_fa_locals = ifelse(
                    p_fa_locals == '1',
                    'none',
                    paste0('\033[0;1m',
                           p_fa_locals,
                           '\033[0m')
                )
                cat(
                    paste0(
                        'Futility bounds: ',
                        paste(
                            paste0('(',
                                   looks[-mlook],
                                   ') ',
                                   p_fa_locals),
                            collapse = '; '
                        ),
                        '\n\033[0;3mLikelihoods of exceeding futility alpha if H0 true: ',
                        paste(paste0(
                            '(',
                            looks[-mlook],
                            ') ',
                            ro(df_pow[[paste0('ratio_futil_', p_nam, '_h0')]][looks[-mlook]], round_to),
                            collapse = '; '
                        )),
                        '\nLikelihoods of exceeding futility alpha if H1 true: ',
                        paste(paste0(
                            '(',
                            looks[-mlook],
                            ') ',
                            ro(df_pow[[paste0('ratio_futil_', p_nam, '_h1')]][looks[-mlook]], round_to),
                            collapse = '; '
                        )),
                        '\033[0m'
                    ),
                    fill = TRUE
                )
            }
        }
    }
    if (length(df_p_names) > 1) {
        toprint = paste0(
            '\033[0;4mGlobal\033[0m ("combined significance") type I error: ',
            ro(df_pow$ratio_combined_sign_h0[df_nrow],
               # same as global_type1
               round_to),
            ' (included: ',
            paste(df_p_names, collapse = ', '),
            '; power for reaching the "combined significance": ',
            ro(df_pow$ratio_combined_sign_h1[df_nrow],
               # same as global_power
               round_to),
            ')'
        )
        if (mlook > 1) {
            toprint = paste0(
                toprint,
                '\n\033[0;3mLikelihoods of stopping for (combined) \033[1msignificance\033[0;3m if H0 true: ',
                paste(paste0(
                    '(',
                    looks,
                    ') ',
                    ro(df_pow$ratio_combined_sign_h0[looks],
                       round_to),
                    collapse = '; '
                )),
                '\nLikelihoods of stopping for (combined) \033[1msignificance\033[0;3m if H1 true: ',
                paste(paste0(
                    '(',
                    looks,
                    ') ',
                    ro(df_pow$ratio_combined_sign_h1[looks],
                       round_to),
                    collapse = '; '
                )),
                '\033[0m'
            )
            if (length(df_fut_names) > 1) {
                toprint = paste0(
                    toprint,
                    '\n\033[0;3mLikelihoods of stopping for (combined) \033[1mfutility\033[0;3m if H0 true: ',
                    paste(paste0(
                        '(',
                        looks[-mlook],
                        ') ',
                        ro(df_pow$ratio_combined_fut_h0[looks[-mlook]],
                           round_to),
                        collapse = '; '
                    )),
                    '\nLikelihoods of stopping for (combined) \033[1mfutility\033[0;3m if H1 true: ',
                    paste(paste0(
                        '(',
                        looks[-mlook],
                        ') ',
                        ro(df_pow$ratio_combined_fut_h1[looks[-mlook]],
                           round_to),
                        collapse = '; '
                    )),
                    '\033[0m'
                )
            }
        }
        cat(toprint, fill = TRUE)
    }
}
