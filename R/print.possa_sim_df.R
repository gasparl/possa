#'@title Print sim results
#'
#'@description Prints information about the simulated p values created by the
#'  \code{\link[POSSA:sim]{POSSA::sim}} function. This is an extension (method)
#'  of the base R \code{\link{print}} function, so it can be called simply as
#'  \code{print()}.
#'@param x The \code{\link{data.frame}} returned by the
#'  \code{\link[POSSA:sim]{POSSA::sim}} function.
#'@param descr_cols When given as a character element or vector, specifies the
#'  factors for which descriptive data should be shown (by group, if
#'  applicable). By default \code{NULL}, it identifies (similar as
#'  \code{group_by}) factors, if any, given to the \code{sim} function (via
#'  \code{fun_obs}) that produced the given \code{x} data.
#'@param descr_func Function used for printing descriptives (see
#'  \code{descr_cols}). By default, it uses the \code{\link{summary}}
#'  (\code{\link{base}}) function.
#'@param group_by When given as a character element or vector, specifies the
#'  factors by which to group the descriptives: the \code{x} data will be
#'  divided into parts by these factors and these parts will be analyzed
#'  separately, with descriptives printed per each part. By default
#'  (\code{NULL}), it identifies factors, if any, given to the \code{sim}
#'  function (via \code{fun_obs}) that created the given \code{x} data.
#'@param ... (Allow additional arguments for technical reasons.)
#'
#'@return Returns nothing (\code{NULL}); this method serves only to print
#'  information to the console.
#'
#'@seealso \code{\link{sim}}
#'
#' @export
print.possa_sim_df = function(x,
                              group_by = NULL,
                              descr_cols = TRUE,
                              descr_func = summary,
                              ...) {
    cat(
        format.possa_sim_df(
            x = x,
            group_by = group_by,
            descr_cols = descr_cols,
            descr_func = descr_func,
            ...
        ),
        fill = TRUE
    )
}

format.possa_sim_df = function(x,
                               group_by = NULL,
                               descr_cols = TRUE,
                               descr_func = summary,
                               ...) {
    .look = NULL
    ._possa_fact_combs = NULL
    cat('\033[0;34mPOSSA sim() results (p values)\033[0m', fill = TRUE)
    df_sim = data.table::copy(x)
    validate_args(match.call(),
                  list(
                      val_arg(df_sim, c('df')),
                      val_arg(group_by, c('null', 'char')),
                      val_arg(descr_cols, c('bool', 'char')),
                      val_arg(descr_func, c('function'), 1)
                  ))
    cat('\033[0;37mSample:\033[0m', fill = TRUE)

    setDT(df_sim)
    print(utils::head(df_sim, min(max(c(
        df_sim$.look, 2L
    )), 8L)))
    setkey(df_sim, .look)

    n_cols = c()
    p_names_all = c()
    fac_cols = c()
    for (c_nam in colnames(df_sim)) {
        col = df_sim[[c_nam]]
        if ('possa_n' %in% class(col)) {
            n_cols = c(n_cols, c_nam)
        } else if ('possa_fac' %in% class(col)) {
            fac_cols = c(fac_cols, c_nam)
        } else if (startsWith(c_nam, 'p_') &&
                   endsWith(c_nam, '_h0')) {
            p_nam_h1 = paste0(substr(c_nam, 1, nchar(c_nam) - 3), '_h1')
            if (p_nam_h1 %in% colnames(df_sim)) {
                p_names_all = c(p_names_all, c_nam)
                p_names_all = c(p_names_all, p_nam_h1)
            }
        }
    }
    if (is.null(group_by)) {
        group_by = fac_cols
    }
    if (length(group_by) > 0) {
        df_sim[, ._possa_fact_combs := do.call(paste, c(.SD, sep = '; ')), .SDcols = group_by]
        possafacts = unique(df_sim$._possa_fact_combs)
        setindex(df_sim, ._possa_fact_combs)
    } else {
        possafacts = NA
    }
    if (descr_cols[1] == TRUE) {
        descr_cols = names(df_sim)[!names(df_sim) %in%
                                       c(
                                           '.iter',
                                           '.look',
                                           '.n_total',
                                           '._possa_fact_combs',
                                           p_names_all,
                                           n_cols,
                                           fac_cols
                                       )]
        if (!length(descr_cols) > 0) {
            descr_cols = FALSE
        }
    }
    if (descr_cols[1] != FALSE) {
        cat('\033[0;37mDescriptives:\033[0m', fill = TRUE)
        for (possa_fact in possafacts) {
            # when none: possa_fact = NA
            if (is.na(possafacts[1])) {
                pvals_df = df_sim
            } else {
                # print descriptives of all included
                if (possa_fact == possafacts[1]) {
                    cat('\033[0;40mTOTAL (all groups)\033[0m',
                        fill = TRUE)
                    for (desc_col in descr_cols) {
                        cat(
                            '\033[0;3m',
                            desc_col,
                            ':\033[0m ',
                            sep = '',
                            fill = TRUE
                        )
                        print(descr_func(df_sim[[desc_col]]))
                    }
                }
                # if applicable, take only given factor combination & print its "group" name
                pvals_df = df_sim[._possa_fact_combs == possa_fact]
                cat(
                    '\nGROUP (',
                    paste(group_by, collapse = '; '),
                    '): \033[0;40m',
                    possa_fact,
                    '\033[0m',
                    fill = TRUE,
                    sep = ''
                )
            }
            for (desc_col in descr_cols) {
                cat('\033[0;3m',
                    desc_col,
                    ':\033[0m ',
                    sep = '',
                    fill = TRUE)
                print(descr_func(pvals_df[[desc_col]]))
            }
        }
    }
}