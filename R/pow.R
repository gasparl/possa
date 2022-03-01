#'@title Power calculation
#'
#'@description Calculates power and local alphas based on simulated p values
#'  (which should be provided as created by the
#'  \code{\link[POSSA:sim]{POSSA::sim}} function). The calculation for
#'  sequential testing involves a staircase procedure during which an initially
#'  provided set of local alphas is continually adjusted until the (approximate)
#'  specified global error rate (e.g., global alpha = .05) is reached: the value
#'  of adjustment is decreasing while global error rate is larger than
#'  specified, and increasing while global error rate is smaller than specified;
#'  a smaller step is chosen whenever the direction (increase vs. decrease)
#'  changes; the procedure stops when the global error rate is close enough to
#'  the specified one (e.g., matches it up to 4 fractional digits) or when the
#'  specified smallest step is passed. The adjustment works via a dedicated
#'  ("\code{adjust}") function that either replaces missing (\code{NA}) values
#'  with varying alternatives or (when there are no missing values) in some
#'  manner varyingly modifies the initial values (e.g. by addition or
#'  multiplication).
#'@param p_values A \code{\link{data.frame}} containing the simulated
#'  iterations, looks, and corresponding H0 and H1 p value outcomes, as returned
#'  by the \code{\link[POSSA:sim]{POSSA::sim}} function. (Custom data frames are
#'  also accepted, but may not work as expected.)
#'@param alpha_locals A number, a numeric vector, or a named \code{\link{list}}
#'  of numeric vectors, that specify the initial set of local alphas that decide
#'  on statistical significance (for interim looks as well as for the final
#'  look), and, if significant, stop the experiment at the given interim look;
#'  to be adjusted via the \code{adjust} function; see the \code{adjust}
#'  parameter below. Any of the numbers included can always be \code{NA} values
#'  as well (which indicates alphas to be calculated; again, see the related
#'  \code{adjust} parameter below). In case of a vector or a list of vectors,
#'  the length of each vector must correspond exactly to the maximum number of
#'  looks in the \code{p_values} data frame. When a \code{\link{list}} is given,
#'  the names of the list element(s) must correspond to the root of the related
#'  H0 and H1 p value column name pair(s) (in the \code{p_values} data frame),
#'  that is, without the "\code{_h0}" and "\code{_h1}" suffixes: for example, if
#'  the column name pair is "\code{p_test4_h0}" and "\code{p_test4_h1}", the
#'  name of the corresponding list element should be "\code{p_test4}". If a
#'  single number or a single numeric vector is given, all potential p value
#'  column pairs are automatically detected as starting with "\code{p_}" prefix
#'  and ending with "\code{_h0}" and "\code{_h1}". In case of a single vector
#'  given, each such automatically detected p value pair receives this same
#'  vector. In case of a single number given, this all elements of all vectors
#'  will be this same number. The default \code{NULL} value specifies "fixed
#'  design" (no interim stopping alphas) with final alpha as specified as
#'  \code{alpha_global}. (This is useful for cases where only futility bounds
#'  are to be set for stopping.)
#'@param alpha_global Global alpha (expected error rate in total); \code{0.05}
#'  by default.
#'@param adjust The function via which the initial vector local alphas is
#'  modified with each step of the staircase procedure. Three arguments are
#'  passed to it: \code{adj}, \code{orig}, and \code{prev}. The \code{adj}
#'  parameter is mandatory; it passes the pivotal changing value that, starting
#'  from an initial value (see \code{adj_init}), should, via the staircase
#'  steps, decrease when the global error rate is too large, and increase when
#'  the global error rate is too small. The \code{orig} parameter (optional)
#'  always passes the same original vector of alphas as they were provided via
#'  \code{alpha_locals}. The \code{prev} parameter (optional) passes the
#'  "latest" vector of local alphas, which were obtained in the previous
#'  adjustment step (or, in the initial run, it is the original vector, i.e.,
#'  the same as \code{orig}). When \code{NULL} (default), a function is given
#'  internally that simply replaces \code{NA}s with the varying adjustment value
#'  (as \code{{ prev[is.na(orig)] = adj; return(prev) }}).
#'@param adj_init The initial adjustment value that is used as the "\code{adj}"
#'  parameter in the "\code{adjust}" function and is continually adjusted via
#'  the staircase steps (see \code{staircase_steps} parameter). When \code{NULL}
#'  (default), it is calculated as the global alpha divided by the maximum
#'  number of looks (Bonferroni correction), as a rough initial approximation
#'  with the assumption that "\code{adj}" is used as a replacement for
#'  \code{NA}s.
#'@param staircase_steps Numeric vector that specifies the (normally decreasing)
#'  sequence of step sizes for the staircase that narrows down on the specified
#'  global error error by decreasing or increasing the adjustment value
#'  (initially: \code{adj_init}): the step size (numeric value) is added for
#'  increase, and subtracted for decrease. Whenever the direction (decrease vs.
#'  increase) is changed, the staircase moves on to the next step size. When the
#'  direction changes and there are no more steps remaining, the procedure is
#'  finished (regardless of the global error rate). By default (\code{NULL}),
#'  the \code{staircase_steps} is either "\code{0.01 * (0.5 ^ (seq(0, 11, 1)))}"
#'  (giving: \code{0.01, 0.005, 0.0025, ...}) or "\code{0.5 * (0.5 ^ (seq(0, 11,
#'  1)))}" (giving: \code{0.05, 0.025, 0.0125, ...}). The latter is chosen when
#'  adjustment via multiplication is assumed, which is simply based on finding
#'  any multiplication sign (\code{*}) in a given custom \code{adjust}
#'  function. The former is chosen in any other case.
#'@param alpha_precision During the staircase procedure, at any point when the
#'  simulated global error rate first matches the given \code{alpha_global} at
#'  least for the number of fractional digits given here
#'  (\code{alpha_precision}; default: \code{5}), the procedure stops and the
#'  results are printed. (Otherwise, the procedures finishes only when all steps
#'  given as \code{staircase_steps} have been used.)
#'@param fut_locals Specifies local futility bounds that may stop the experiment
#'  at the given interim looks if the corresponding p value is above the given
#'  futility bound value. When \code{NULL} (default), sets no futility bounds.
#'  Otherwise, it follows the same logic as \code{alpha_locals} and has the same
#'  input possibilities (number, numeric vector, or named list of numeric
#'  vectors).
#'@param multi_logic When multiple p values are evaluated for stopping rules,
#'  \code{multi_logic} specifies the function used for how to evaluate the
#'  multiple outcomes as a single \code{TRUE} or \code{FALSE} value that decides
#'  whether or not to stop at a given look. The default, \code{'all'}, specifies
#'  that all of the p values must pass the boundary for stopping. The other
#'  acceptable character input is \code{'any'}, which specifies that the
#'  collection stops when any of the p values pass the boundary for stopping.
#'  Instead of these strings, the actual \code{\link{all}} and \code{\link{any}}
#'  would lead to identical outcomes, respectively, but the processing would be
#'  far slower (since the string \code{'all'} or \code{'any'} inputs specify a
#'  dedicated faster internal solution). For custom combinations, any custom
#'  function can be given, which will take, as arguments, the p value columns in
#'  their given order (either in the \code{p_values} data frame, or as specified
#'  in \code{alpha_locals}), and should return a single \code{TRUE} or
#'  \code{FALSE} value.
#'@param multi_logic_fut Same as \code{multi_logic}, but for futility bounds
#'  (for the columns specified in \code{fut_locals}).
#'@param group_by When given as a character element or vector, specifies the
#'  factors by which to group the analysis: the \code{p_values} data will be
#'  divided into parts by these factors and these parts will be analyzed
#'  separately, with power and error information (and descriptives, if
#'  specified) printed per each part. By default (\code{NULL}), it identifies
#'  factors, if any, given to the \code{sim} function (via \code{fun_obs}) that
#'  created the given \code{p_values} data.
#'@param alpha_locals_extra Optional extra and "non-stopper" alphas via which to
#'  evaluate p values per look, but without stopping the data collection
#'  regardless of statistical significance.
#'@param design_fix Whether to calculate fixed design (as comparable
#'  alternative(s) to given "looks" of the sequential design). If \code{NULL}
#'  (default), shows power for maximum sample (last "look") only. If set to
#'  \code{TRUE}, shows power for each given "look" (interim and final).
#'  Altogether omitted when set to \code{FALSE}.
#'@param design_seq Whether to calculate sequential design (default:
#'  \code{TRUE}). (Although this R package is designed for sequential analysis,
#'  it can be used for fixed designs alone too.) Altogether omitted when set to
#'  \code{FALSE}.
#'@param descr_cols When given as a character element or vector, specifies the
#'  factors for which descriptive data should be shown (by group, if
#'  applicable). By default (\code{TRUE}), it identifies (similar as
#'  \code{group_by}) factors, if any, given to the \code{sim} function (via
#'  \code{fun_obs}) that produced the given \code{p_values} data. If set to
#'  \code{FALSE}, no descriptive data is shown.
#'@param descr_func Function used for printing descriptives (see
#'  \code{descr_cols}). By default, it uses the \code{\link{summary}}
#'(\code{\link{base}}) function.
#'@param round_to Number of fractional digits (default: \code{5}) to round to,
#'  for the displayed power and error rate information.
#'@param seed Number for \code{\link{set.seed}}; \code{8} by default. Set to
#'  \code{NULL} for random seed.
#'
#'@return Apart from printing the main power, error rate, and adjusted local
#'  alpha outcomes, the function (invisibly) returns a \code{\link{data.frame}}
#'  that includes all details of the calculated information.
#'
#'@note
#'
#'@note
#'
#'For the replicability, in case the \code{adjust} function uses any
#'randomization, \code{\link{set.seed}} is executed in the beginning of this
#'function, each time it is called; see the \code{seed} parameter.
#'
#'This function uses, internally, the \code{\link{data.table}} R package.
#'
#'@references
#'
#'Lakens, D. (2014). Performing high-powered studies efficiently with sequential
#'analyses: Sequential analyses. European Journal of Social Psychology, 44(7),
#'701–710. \doi{https://doi.org/10.1002/ejsp.2023}
#'
#'@seealso \code{\link{sim}}
#' @examples
#' # some pow
#'
#' @export
pow = function(p_values,
               alpha_locals = NULL,
               alpha_global = 0.05,
               adjust = NULL,
               adj_init = NULL,
               staircase_steps = NULL,
               alpha_precision = 5,
               fut_locals = NULL,
               multi_logic = 'all',
               multi_logic_fut = 'all',
               group_by = NULL,
               alpha_locals_extra = NULL,
               design_fix = NULL,
               design_seq = TRUE,
               descr_cols = TRUE,
               descr_func = summary,
               round_to = 5,
               seed = 8) {
    validate_args(
        match.call(),
        list(
            val_arg(p_values, c('df')),
            val_arg(alpha_global, c('num'), 1),
            val_arg(adjust, c('null', 'function')),
            val_arg(adj_init, c('null', 'num'), 1),
            val_arg(group_by, c('null', 'char')),
            val_arg(design_fix, c('null', 'bool'), 1),
            val_arg(design_seq, c('bool'), 1),
            val_arg(descr_cols, c('bool'), 1),
            val_arg(descr_func, c('function'), 1),
            val_arg(round_to, c('num'), 1),
            val_arg(multi_logic, c('char', 'num'), 1, c('all', 'any')),
            val_arg(multi_logic_fut, c('char', 'num'), 1, c('all', 'any')),
            val_arg(staircase_steps, c('null', 'num')),
            val_arg(alpha_precision, c('num'), 1)
        )
    )
    look = NULL
    iter = NULL
    h0_stoP = NULL
    h1_stoP = NULL
    h0_stoP_fa = NULL
    h1_stoP_fa = NULL
    min_look = NULL
    ._possa_fact_combs = NULL
    . = NULL

    set.seed(seed)
    if (!'possa_df' %in% class(p_values)) {
        warning(
            'The given data frame seems not to have been created by the ',
            '"possa::sim_pvals()" function; it may not fit the "possa::get_pow()" function.',
            immediate. = TRUE
        )
    }
    if (is.function(multi_logic)) {
        if (isTRUE(all.equal(multi_logic, any))) {
            message(
                'Warning: Operation is hugely faster with ',
                '"any" string argument for multi_logic.'
            )
        } else if (isTRUE(all.equal(multi_logic, all))) {
            message(
                'Warning: Operation is hugely faster with ',
                '"all" string argument for multi_logic.'
            )
        }
    } else {
        if (multi_logic == 'all') {
            multi_logic = `|`
        } else if (multi_logic == 'any') {
            multi_logic = `&`
        }
    }
    if (is.function(multi_logic_fut)) {
        if (isTRUE(all.equal(multi_logic_fut, any))) {
            message(
                'Warning: Operation is hugely faster with ',
                '"any" string argument for multi_logic_fut.'
            )
        } else if (isTRUE(all.equal(multi_logic_fut, all))) {
            message(
                'Warning: Operation is hugely faster with ',
                '"all" string argument for multi_logic_fut.'
            )
        }
    } else {
        if (multi_logic_fut == 'all') {
            multi_logic_fut = `|`
        } else if (multi_logic_fut == 'any') {
            multi_logic_fut = `&`
        }
    }
    setDT(p_values)
    setkey(p_values, look)
    setindex(p_values, iter)
    looks = unique(p_values$look)
    mlook = max(p_values$look)
    # get columns with sample sizes (n), factors (fac), and p values (p_/_h0/1)
    n_cols = c()
    fac_cols = c()
    p_names_auto = c()
    for (c_nam in colnames(p_values)) {
        col = p_values[[c_nam]]
        if ('possa_n' %in% class(col)) {
            n_cols = c(n_cols, c_nam)
        } else if ('possa_fac' %in% class(col)) {
            fac_cols = c(fac_cols, c_nam)
        } else if (startsWith(c_nam, 'p_') &&
                   endsWith(c_nam, '_h0')) {
            pnam = substr(c_nam, 1, nchar(c_nam) - 3)
            if (paste0(pnam, '_h1') %in% colnames(p_values)) {
                p_names_auto = c(p_names_auto, pnam)
            }
        }
    }
    a_locals = list()
    p_names = p_names_auto
    # extract (if given) predetermined local alphas and/or specified p value columns
    if (!is.null(alpha_locals)) {
        if (is.atomic(alpha_locals)) {
            # if vector given
            if (is.character(alpha_locals)) {
                # if character, simply give the names
                loc_pnames = alpha_locals
            } else {
                # if numeric, assign to each p column
                if (length(alpha_locals) == 1) {
                    for (pnam in p_names) {
                        a_locals[[pnam]] = rep(alpha_locals, mlook)
                    }
                } else if (!length(alpha_locals) == mlook) {
                    stop(
                        'Wrong argument for "alpha_locals". (If a numeric vector is given, ',
                        'it must have same length as the maximum number of looks (in this case ',
                        mlook,
                        ').)'
                    )
                } else {
                    for (pnam in p_names) {
                        a_locals[[pnam]] = alpha_locals
                    }
                }
                loc_pnames = p_names
            }
        } else {
            # if list, use as it is, and assign per name
            for (a_vec in alpha_locals) {
                if (!(is.atomic(a_vec) && length(a_vec) == mlook)) {
                    stop(
                        'Wrong argument for "alpha_locals". (If a list is given, ',
                        'it must consist of one or more vectors of numbers with the',
                        ' same length as the maximum number of looks (in this case ',
                        mlook,
                        ').)'
                    )
                }
            }
            a_locals = alpha_locals
            loc_pnames = names(alpha_locals)
        }
        for (pname in loc_pnames) {
            if (!pname %in% p_names) {
                stop(
                    'Wrong argument for "alpha_locals". ',
                    'There is no column name pair "',
                    pname,
                    '0"/"',
                    pname,
                    '1".'
                )
            }
        }
        p_names = loc_pnames
    }
    if (!length(a_locals) > 0) {
        # if not given, check only last look (with alpha_global)
        # for the rest, local alpha will be 0
        for (pnam in p_names) {
            a_locals[[pnam]] = c(rep(0, (mlook - 1)), alpha_global)
        }
    } else {
        lapply(a_locals, function(vec) {
            vec = vec[!is.na(vec)]
            if (any(vec > 1 | vec < 0)) {
                stop(
                    'All alpha values given in alpha_locals must be',
                    ' between 0 and 1 (or NA).'
                )
            }
        })
    }
    p_extr = NULL
    if (!is.null(alpha_locals_extra)) {
        lapply(alpha_locals_extra, function(vec) {
            if (any(is.na(vec))) {
                stop(
                    'The alpha values given in alpha_locals_extra must ',
                    'not contain NA values.'
                )
            } else if (any(vec > 1 | vec < 0)) {
                stop(
                    'All alpha values given in alpha_locals_extra must be',
                    ' between 0 and 1.'
                )
            } else if (!(is.atomic(a_vec) &&
                         length(a_vec) == mlook)) {
                stop(
                    'Wrong argument for "alpha_locals_extra". (If a list is given, ',
                    'it must consist of one or more vectors of numbers with the',
                    ' same length as the maximum number of looks (in this case ',
                    mlook,
                    ').)'
                )
            }
        })
        p_extr = names(alpha_locals_extra)
        for (pname in p_extr) {
            if (!pname %in% p_names_auto) {
                stop(
                    'Wrong argument for "alpha_locals_extra". ',
                    'There is no column name pair "',
                    pname,
                    '0"/"',
                    pname,
                    '1".'
                )
            }
        }
        a_locals = append(a_locals, alpha_locals_extra)
    }
    p_names_extr = c(p_names, p_extr)
    fa_locals = list()
    # extract (if given) predetermined futility bounds
    if (!is.null(fut_locals)) {
        if (is.atomic(fut_locals)) {
            # if vector given, assign to each p column
            if (length(fut_locals) == 1) {
                for (pnam in p_names) {
                    fa_locals[[pnam]] = rep(fut_locals, mlook)
                }
            } else if (!length(fut_locals) == (mlook - 1)) {
                stop(
                    'Wrong argument for "fut_locals". (If a numeric vector is given, ',
                    'its length must be one less than the maximum number of looks (in this case ',
                    (mlook - 1),
                    ').)'
                )
            } else {
                for (pnam in p_names) {
                    fa_locals[[pnam]] = fut_locals
                }
            }
        } else {
            # if list, assign per name
            for (a_vec in fut_locals) {
                if (!(is.atomic(a_vec) && length(a_vec) == (mlook - 1))) {
                    stop(
                        'Wrong argument for "fut_locals". (If a list is given, ',
                        'it must consist of one or more vectors of numbers with a ',
                        'length one less than the maximum number of looks (in this case ',
                        mlook,
                        ').)'
                    )
                }
                fa_locals = fut_locals
            }
            for (pname in names(fut_locals)) {
                if (!pname %in% p_names) {
                    stop(
                        'Wrong argument for "fut_locals". ',
                        'There is no column name pair "',
                        pname,
                        '0"/"',
                        pname,
                        '1".'
                    )
                }
            }
        }
    }
    # if not given, add 1 for all local futility bounds
    if (!length(fa_locals) > 0) {
        for (pnam in p_names) {
            fa_locals[[pnam]] = rep(1, (mlook - 1))
        }
    }


    if (is.null(staircase_steps)) {
        if (anyNA(unlist(a_locals))) {
            # default for NA replacement: 11 steps from 0.01, decreasing by halves
            # check: formatC(staircase_steps, digits = 12, format = "f")
            staircase_steps = 0.01 * (0.5 ** (seq(0, 11, 1)))
        } else if (!is.null(adjust)) {
            staircase_steps = 0.01 * (0.5 ** (seq(0, 11, 1)))
            if (any(grepl('*', deparse(body(adjust)), fixed = TRUE))) {
                if (is.null(adj_init)) {
                    adj_init = 1 # assumes multiplication
                }
                # slightly larger steps, again assuming multiplications
                staircase_steps = 0.5 * (0.5 ** (seq(0, 11, 1)))
            }
        } else {
            staircase_steps = NA
        }
    } else if (any(is.na(staircase_steps))) {
        stop('The staircase_steps vector must not contain NA values.')
    }

    if (is.null(adj_init)) {
        # if not given, start adjustment with bonferroni correction
        # (assuming it's going to be used as an alpha level, replacing NAs)
        adj_init = alpha_global / mlook
    }
    if (is.null(adjust)) {
        adjust = function(adj, prev, orig) {
            prev[is.na(orig)] = adj
            return(prev)
        }
    } else {
        if (!'adj' %in% methods::formalArgs(adjust)) {
            stop('The "adjust" function must contain an "adj" parameter.')
        }
        for (a_arg in c('prev', 'orig')) {
            if (!a_arg %in%  methods::formalArgs(adjust)) {
                formals(adjust)[[a_arg]] = NA
            }
        }
        a_example = a_locals[p_names[1]]
        adjusted_check = adjust(adj = adj_init,
                                prev = a_example,
                                orig = a_example)
        if (length(adjusted_check) != mlook) {
            stop(
                'The "adjust" function must return a vector with the same length ',
                ' as the maximum number of looks (in this case ',
                mlook,
                ').'
            )
        }
        if (anyNA(adjusted_check)) {
            warning(
                'The local alphas returned by the given "adjust" contain NA(s). ',
                'This should normally not happen, and may cause errors.'
            )
        }
        if (sum(adjusted_check, na.rm = TRUE) > sum(adjust(
            adj = adj_init + 0.1,
            prev = a_example,
            orig = a_example
        ),
        na.rm = TRUE) |
        sum(adjusted_check, na.rm = TRUE) < sum(adjust(
            adj = adj_init - 0.1,
            prev = a_example,
            orig = a_example
        ),
        na.rm = TRUE)) {
            warning(
                'The local alphas returned by the  "adjust" function should normally ',
                'increase when the "adj" argument increases, and decrease when the latter ',
                'decreases. The behavior of the given function seems different.'
            )
        }
    }


    if (is.null(group_by)) {
        group_by = fac_cols
    } else if (!identical(sort(group_by), sort(fac_cols))) {
        message('Custom "group_by" argument given. Be cautious.')
    }
    if (length(group_by) > 0) {
        p_values[, ._possa_fact_combs := do.call(paste, c(.SD, sep = '; ')), .SDcols = group_by]
        possafacts = unique(p_values$._possa_fact_combs)
    } else {
        possafacts = NA
    }
    p_names_h0 = paste0(p_names, '_h0')
    p_names_h1 = paste0(p_names, '_h1')
    if (descr_cols[1] == TRUE) {
        descr_cols = names(p_values)[!names(p_values) %in%
                                         c(
                                             'iter',
                                             'look',
                                             paste0(p_names_extr, '_h0'),
                                             paste0(p_names_extr, '_h1'),
                                             n_cols,
                                             fac_cols
                                         )]
        if (!length(descr_cols) > 0) {
            descr_cols = FALSE
        }
    }
    out_dfs = list()
    # calculate results separately for each factor combination
    for (possa_fact in possafacts) {
        # when none: possa_fact = NA
        if (is.na(possafacts)) {
            pvals_df = p_values
        } else {
            # print descriptives of all included
            if (descr_cols[1] != FALSE &&
                possa_fact == possafacts[1]) {
                cat('-- DESCRIPTIVES (total) --', fill = TRUE)
                for (desc_col in descr_cols) {
                    cat(desc_col, ': ', sep = '')
                    descr_func(p_values[[desc_col]])
                }
            }
            # if applicable, take only given factor combination & print its "group" name
            pvals_df = p_values[._possa_fact_combs == possa_fact]
            cat('GROUP: ', possa_fact, fill = TRUE)
        }
        if (descr_cols[1] != FALSE) {
            cat('-- DESCRIPTIVES --', fill = TRUE)
            for (desc_col in descr_cols) {
                cat(desc_col,
                    ': ',
                    sep = '',
                    fill = TRUE)
                print(descr_func(pvals_df[[desc_col]]))
            }
            cat('', fill = TRUE)
        }
        tot_samples = c()
        for (lk in looks) {
            tot_samples = c(tot_samples, sum(pvals_df[.(lk), .SD, .SDcols = n_cols, mult = 'first']))
        }
        look_ratios = tot_samples / tot_samples[mlook]
        ## Fixed design calculation below
        if (is.null(design_fix)) {
            fix_looks = mlook # (default) show at max look only
        } else  if (design_fix == TRUE) {
            fix_looks = 1:mlook # show outcome at all looks
        } else if (design_fix == FALSE) {
            fix_looks = NULL # show none
        }
        for (f_look in fix_looks) {
            pvals_df_fix = pvals_df[look == f_look]
            cat('-- FIXED DESIGN; N(total) =',
                tot_samples[f_look],
                '--',
                fill = TRUE)
            for (p_nam in p_names_extr) {
                cat(
                    '(',
                    p_nam,
                    ') Type I error: ',
                    round(mean(
                        pvals_df_fix[[paste0(p_nam, '_h0')]] < alpha_global
                    ), round_to),
                    '; Power: ',
                    round(mean(
                        pvals_df_fix[[paste0(p_nam, '_h1')]] < alpha_global
                    ), round_to),
                    '\n',
                    sep = '',
                    fill = TRUE
                )
            }
        }

        ## Sequential design calculation below (when applicable)
        if (design_seq == TRUE & mlook > 1) {
            # "trial and error" straircase procedure below, to get the desired adjusted alpha
            locls_temp = a_locals
            a_adj = adj_init
            stair_steps = staircase_steps
            if (!is.null(alpha_locals_extra) &
                !is.na(stair_steps[length(stair_steps)])) {
                stair_steps = c(stair_steps, NA)
            }
            a_step = stair_steps[1] # initial a_adj-changing step
            p_h0_sign_names = paste0(p_names, '_h0_sign')
            multi_p = length(p_h0_sign_names) > 1
            p_h0_sign_names_plus = c(p_h0_sign_names, 'look', 'iter')
            p_h0_fut_names = paste0(p_names, '_h0_fut')
            p_h1_sign_names = paste0(p_names, '_h1_sign')
            p_h1_sign_names_plus = c(p_h1_sign_names, 'look', 'iter')
            p_h1_fut_names = paste0(p_names, '_h1_fut')
            pb = utils::txtProgressBar(
                min = 0,
                max = length(staircase_steps),
                initial = 0,
                style = 3
            )
            p_names_temp = p_names
            for (p_nam in p_names_temp) {
                pvals_df[, c(paste0(p_nam, '_h0_sign')) := NA] # create sign column for given p
                if (!is.null(fut_locals)) {
                    # if futility bounds are given
                    pvals_df[, c(paste0(p_nam, '_h0_fut')) := TRUE] # create fut column for given p
                }
            }
            safe_count = 0
            type1 = 1
            while (length(stair_steps) > 0) {
                ### TO REMOVE (just for testing)
                # cat('\ntype1',
                #     type1,
                #     'a_adj:',
                #     a_adj,
                #     'new step:',
                #     a_step,
                #     fill = T)

                # calculate H0 significances (T/F) & stops (T/F) based on adjusted alphas
                if (is.na(stair_steps[1])) {
                    # as a last step, add non-stopping columns, if any
                    p_names_temp = p_names_extr
                }
                safe_count = safe_count + 1
                if (safe_count %% 100 == 0) {
                    cat(
                        paste0(
                            'There have been ',
                            safe_count,
                            ' iterations. The arguments given may be problematic.',
                            '\nNote: The current global type 1 error is ',
                            as.numeric(round(type1, 4)),
                            ' (expected: ',
                            round(alpha_global, 4),
                            ').'
                        ),
                        fill = TRUE
                    )
                    if (readline('Do you want to quit this function? (y/n) \n') == 'y') {
                        return(cat("Exited safely."))
                    }
                }
                for (p_nam in p_names_temp) {
                    # adjust alpha
                    locls_temp[[p_nam]] = adjust(adj = a_adj,
                                                 prev = locls_temp[[p_nam]],
                                                 orig = a_locals[[p_nam]])
                    for (lk in 1:mlook) {
                        # decide significance at given look for given p
                        pvals_df[.(lk), c(paste0(p_nam, '_h0_sign')) :=
                                     .SD < locls_temp[[p_nam]][lk], .SDcol = paste0(p_nam, '_h0')]
                        if (!is.null(fut_locals) &
                            lk != mlook) {
                            # (futility never matters at max look)
                            pvals_df[.(lk), c(paste0(p_nam, '_h0_fut')) :=
                                         .SD > fa_locals[[p_nam]][lk], .SDcol = paste0(p_nam, '_h0')]
                        }
                    }
                }

                # now check the global type 1 error

                # first check at which look we stop
                # if multiple p values are given as stoppers
                # use multi_logic to check where to stop
                # otherwise it simply stops where the given p is sign
                if (multi_p) {
                    if (is.function(multi_logic)) {
                        pvals_df[,  h0_stoP := apply(.SD, 1, multi_logic), .SDcols = p_h0_sign_names]
                    } else {
                        pvals_df[, h0_stoP := Reduce(multi_logic, .SD), .SDcols = p_h0_sign_names]
                    }
                } else {
                    pvals_df[, h0_stoP := .SD, .SDcols = p_h0_sign_names]
                }
                if (!is.null(fut_locals)) {
                    # analogue with futility bounds
                    if (multi_p) {
                        if (is.function(multi_logic_fut)) {
                            pvals_df[,  h0_stoP_fa := apply(.SD, 1, multi_logic_fut), .SDcols = p_h0_fut_names]
                        } else {
                            pvals_df[, h0_stoP_fa := Reduce(multi_logic_fut, .SD), .SDcols = p_h0_fut_names]
                        }
                    } else {
                        pvals_df[, h0_stoP_fa := .SD, .SDcols = p_h0_fut_names]
                    }
                    pvals_df[, h0_stoP := h0_stoP | h0_stoP_fa]
                }

                # now get all outcomes at stopping point
                pvals_stp = pvals_df[look == mlook |
                                         h0_stoP == TRUE,  .SD, .SDcols = p_h0_sign_names_plus]
                type1 = mean(unlist(pvals_stp[, min_look := min(look), by = iter][look == min_look, .SD, .SDcols = p_h0_sign_names]))

                if (is.na(stair_steps[1])) {
                    # NA indicates no stairs; nothing left to be done
                    break
                } else{
                    # break if global alpha is correct
                    # otherwise continue the staircase
                    if (round(type1, alpha_precision) == round(alpha_global, alpha_precision)) {
                        if (is.null(alpha_locals_extra)) {
                            break
                        } else {
                            stair_steps = NA
                            a_step = 0
                        }
                    } else if ((type1 < alpha_global &&
                                a_step < 0) ||
                               (type1 > alpha_global &&
                                a_step > 0)) {
                        # change staircase direction (and also decrease step) if needed
                        a_step = -stair_steps[1] * sign(a_step)
                        stair_steps = stair_steps[-1]
                        utils::setTxtProgressBar(pb, utils::getTxtProgressBar(pb) + 1)
                    }
                    # adjust a_adj itself by the step
                    # this a_adj will be used in the adjust function
                    a_adj = a_adj + a_step
                }
            }
            utils::setTxtProgressBar(pb, length(staircase_steps))
            close(pb)
            # assign final local alphas
            a_locals_fin = locls_temp

            # calculate H1 significances (T/F) & stops (T/F) based on final alphas
            # (analogue of the H0 significances above)
            for (p_nam in p_names_extr) {
                pvals_df[, c(paste0(p_nam, '_h1_sign')) := NA]
                if (!is.null(fut_locals)) {
                    pvals_df[, c(paste0(p_nam, '_h1_fut')) := TRUE]
                }
                for (lk in 1:mlook) {
                    pvals_df[.(lk), c(paste0(p_nam, '_h1_sign')) :=
                                 .SD < locls_temp[[p_nam]][lk], .SDcol = paste0(p_nam, '_h1')]
                    if (!is.null(fut_locals) &
                        lk != mlook) {
                        # (futility never matters at max look)
                        pvals_df[.(lk), c(paste0(p_nam, '_h1_fut')) :=
                                     .SD > fa_locals[[p_nam]][lk], .SDcol = paste0(p_nam, '_h1')]
                    }
                }
            }
            # now check the global power
            # if multiple p columns, check at which look we stop
            if (multi_p) {
                if (is.function(multi_logic)) {
                    pvals_df[, h1_stoP := apply(.SD, 1, multi_logic), .SDcols = p_h1_sign_names]
                } else {
                    pvals_df[, h1_stoP := Reduce(multi_logic, .SD), .SDcols = p_h1_sign_names]
                }
            } else {
                pvals_df[, h1_stoP := .SD, .SDcols = p_h1_sign_names]
            }

            if (!is.null(fut_locals)) {
                if (multi_p) {
                    if (is.function(multi_logic_fut)) {
                        pvals_df[,  h1_stoP_fa := apply(.SD, 1, multi_logic_fut), .SDcols = p_h1_fut_names]
                    } else {
                        pvals_df[, h1_stoP_fa := Reduce(multi_logic_fut, .SD), .SDcols = p_h1_fut_names]
                    }
                } else {
                    pvals_df[, h1_stoP_fa := .SD, .SDcols = p_h1_fut_names]
                }
                pvals_df[, h1_stoP := h1_stoP | h1_stoP_fa]
            }
            if (multi_p) {
                # if multi_p, get global power at stopping point
                pvals_stp = pvals_df[look == mlook |
                                         h1_stoP == TRUE,  .SD, .SDcols = p_h1_sign_names_plus]
                seq_power = mean(unlist(pvals_stp[, min_look := min(look), by = iter][look == min_look, .SD, .SDcols = p_h1_sign_names]))
            }

            # calculate sample size information per look
            ps_sub0 = data.table::copy(pvals_df)
            ps_sub1 = data.table::copy(pvals_df)
            iters_tot = length(unique(pvals_df$iter))
            stops = list() # collect info per each stop
            previous_h0 = iters_tot # start with max for both
            previous_h1 = iters_tot
            for (lk in looks) {
                # get iterations stopped at given look
                iters_out0 = ps_sub0[look == lk &
                                         h0_stoP == TRUE]
                # remove stopped iterations
                ps_sub0 = ps_sub0[!iter %in% iters_out0$iter,]
                # (same for H1)
                iters_out1 = ps_sub1[look == lk &
                                         h1_stoP == TRUE, ]
                ps_sub1 = ps_sub1[!iter %in% iters_out1$iter, ]
                outs = c()
                # get info per p value column
                for (p_nam in p_names_extr) {
                    # number of significant findings
                    outs[paste0('iters_sign_', p_nam, '_h0')] =
                        sum(iters_out0[[paste0(p_nam, '_h0_sign')]])
                    outs[paste0('iters_sign_', p_nam, '_h1')] =
                        sum(iters_out1[[paste0(p_nam, '_h1_sign')]])
                    if (!is.null(fut_locals)) {
                        # number of futility bound crossings
                        outs[paste0('iters_futil_', p_nam, '_h0')] =
                            sum(iters_out1[[paste0(p_nam, '_h0_fut')]])
                        outs[paste0('iters_futil_', p_nam, '_h1')] =
                            sum(iters_out1[[paste0(p_nam, '_h1_fut')]])
                    }
                }
                # whatever remained
                outs['iters_remain_h0'] =
                    length(unique(ps_sub0$iter))
                outs['iters_remain_h1'] =
                    length(unique(ps_sub1$iter))
                if (lk == mlook) {
                    # at last look, all is stopped that previously remained
                    outs['iters_stopped_h0'] = previous_h0
                    outs['iters_stopped_h1'] = previous_h1
                } else {
                    # calculate stops as previous minus current remaining
                    outs['iters_stopped_h0'] = previous_h0 - outs['iters_remain_h0']
                    outs['iters_stopped_h1'] = previous_h1 - outs['iters_remain_h1']
                }
                stops[[length(stops) + 1]] = c(look = lk,
                                               sample = tot_samples[lk],
                                               outs)
                # assign current remaining as the next "previous remaining"
                previous_h0 = outs['iters_remain_h0']
                previous_h1 = outs['iters_remain_h1']
            }
            df_stops = as.data.frame(do.call(rbind, stops))
            # derive end info (type 1, power, etc) from the iteration ratios
            for (pnam in p_names_extr) {
                # write out local alphas
                df_stops[paste0('alpha_local_', p_nam)] = a_locals_fin[[pnam]]
                if (!is.null(fut_locals)) {
                    # write out local futility bounds (if any)
                    df_stops[paste0('futil_local_', p_nam)] = c(fa_locals[[pnam]], NA)
                }
                for (h01 in c('_h0', '_h1')) {
                    # ratio of significant findings at any given stop
                    df_stops[paste0('ratio_sign_', p_nam, h01)] = df_stops[paste0('iters_sign_', p_nam, h01)] / iters_tot
                    if (!is.null(fut_locals)) {
                        # ratio of futility bound crossings at any given stop
                        df_stops[paste0('ratio_futil_', p_nam, h01)] = df_stops[paste0('iters_futil_', p_nam, h01)] / iters_tot
                    }
                }
            }
            # count ratios of stops per each look
            df_stops$ratio_stopped_h0 = df_stops$iters_stopped_h0 / iters_tot
            df_stops$ratio_stopped_h1 = df_stops$iters_stopped_h1 / iters_tot
            # count cumulative ratios of remainings
            df_stops$ratio_remain_h0 = df_stops$iters_remain_h0 / iters_tot
            df_stops$ratio_remain_h1 = df_stops$iters_remain_h1 / iters_tot
            # the average sample proportion at the given look (mainly just to calculate average-total)
            df_stops$samp_avg_prop_0 = df_stops$ratio_stopped_h0 * df_stops$sample
            df_stops$samp_avg_prop_1 = df_stops$ratio_stopped_h1 * df_stops$sample
            # calculate sums per each column (with different meaning in each case)
            df_stops = rbind(df_stops, colSums(df_stops))
            df_nrow = nrow(df_stops)
            df_stops$look[df_nrow] = 'totals'

            # print results for sequential design
            cat(
                '-- SEQUENTIAL DESIGN; N(average-total) = ',
                round(df_stops$samp_avg_prop_0[df_nrow], 2),
                ' (if H0 true) or ',
                round(df_stops$samp_avg_prop_1[df_nrow], 2),
                ' (if H1 true) --',
                sep = '',
                fill = TRUE
            )

            fut_text = ''
            for (p_nam in p_names_extr) {
                if (!is.null(fut_locals)) {
                    fut_text = paste('\nFutility bounds:',
                                     paste(
                                         paste0(
                                             '(',
                                             looks[-mlook],
                                             ') ',
                                             round(fa_locals[[pnam]], round_to)
                                         ),
                                         collapse = '; '
                                     ))
                }
                cat(
                    '(',
                    p_nam,
                    ') Type I error: ',
                    round(df_stops[[paste0('ratio_sign_', p_nam, '_h0')]][df_nrow], round_to),
                    '; Power: ',
                    round(df_stops[[paste0('ratio_sign_', p_nam, '_h1')]][df_nrow], round_to),
                    '\nAdjusted local alphas: ',
                    paste(paste0(
                        '(',
                        looks,
                        ') ',
                        round(a_locals_fin[[pnam]], round_to)
                    ), collapse = '; '),
                    fut_text,
                    sep = '',
                    fill = TRUE
                )
            }
            if (multi_p) {
                cat(
                    'Global (average) type I error: ',
                    round(type1, round_to),
                    ' (included: ',
                    paste(p_names, collapse = ', '),
                    '; average power: ',
                    round(seq_power, round_to),
                    ')',
                    sep = '',
                    fill = TRUE
                )
            }
            out_dfs[[length(out_dfs) + 1]] = df_stops
        }
    }
    if (design_seq == TRUE & mlook > 1) {
        if (is.na(possafacts)) {
            out_dfs = out_dfs[[1]]
        }
        invisible(out_dfs)
    } else {
        invisible(NULL)
    }
}
