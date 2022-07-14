#'@title Power calculation
#'
#'@description Calculates power and local alphas based on simulated p values
#'  (which should be provided as created by the
#'  \code{\link[POSSA:sim]{POSSA::sim}} function). The calculation for
#'  sequential testing involves a staircase procedure during which an initially
#'  provided set of local alphas is continually adjusted until the (approximate)
#'  specified global type 1 error rate (e.g., global alpha = .05) is reached:
#'  the value of adjustment is decreasing while global type 1 error rate is
#'  larger than specified, and increasing while global type 1 error rate is
#'  smaller than specified; a smaller step is chosen whenever the direction
#'  (increase vs. decrease) changes; the procedure stops when the global type 1
#'  error rate is close enough to the specified one (e.g., matches it up to 4
#'  fractional digits) or when the specified smallest step is passed. The
#'  adjustment works via a dedicated ("\code{adjust}") function that either
#'  replaces missing (\code{NA}) values with varying alternatives or (when there
#'  are no missing values) in some manner varyingly modifies the initial values
#'  (e.g. by addition or multiplication).
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
#'  vector. In case of a single number given, all elements of all vectors will
#'  be assigned this same number (up to the maximum number of looks). If a list
#'  is given and any of the elements contain just one number, it will be
#'  extended into a vector (up to the maximum number of looks). The default
#'  \code{NULL} value specifies "fixed design" (no interim stopping alphas) with
#'  final alpha as specified as \code{alpha_global}, without adjustment
#'  procedure as long as the \code{adjust} argument is also left as default
#'  \code{TRUE}. (This is useful for cases where only futility bounds are to be
#'  set for stopping.)
#'@param alpha_global Global alpha (expected type 1 error rate in total);
#'  \code{0.05} by default. See also \code{multi_logic_global} for when multiple
#'  p values are being evaluated.
#'@param adjust The function via which the initial vector local alphas is
#'  modified with each step of the staircase procedure. Three arguments are
#'  passed to it: \code{adj}, \code{orig}, and \code{prev}. The \code{adj}
#'  parameter is mandatory; it passes the pivotal changing value that, starting
#'  from an initial value (see \code{adj_init}), should, via the staircase
#'  steps, decrease when the global type 1 error rate is too large, and increase
#'  when the global type 1 error rate is too small. The \code{orig} parameter
#'  (optional) always passes the same original vector of alphas as they were
#'  provided via \code{alpha_locals}. The \code{prev} parameter (optional)
#'  passes the "latest" vector of local alphas, which were obtained in the
#'  previous adjustment step (or, in the initial run, it is the original vector,
#'  i.e., the same as \code{orig}). When \code{TRUE} (default), if the given
#'  \code{alpha_locals} contains any \code{NA}s, an \code{adjust} function is
#'  given internally that simply replaces \code{NA}s with the varying adjustment
#'  value (as \code{{ prev[is.na(orig)] = adj; return(prev) }}). If
#'  \code{alpha_locals} contains no \code{NA}s, an \code{adjust} function is
#'  given that multiplies each original local alpha with the varying adjustment
#'  value (as \code{{ return(orig * adj) }}). When set to \code{FALSE}, there
#'  will be no adjustment (staircase procedure omitted): this is useful to
#'  calculate the global type 1 error rate for any given set of local alphas.
#'  Furthermore, if both \code{adjust} and \code{alpha_locals} are left as
#'  default (\code{TRUE} and \code{NULL}), the staircase procedure will be
#'  omitted.
#'@param adj_init The initial adjustment value that is used as the "\code{adj}"
#'  parameter in the "\code{adjust}" function and is continually adjusted via
#'  the staircase steps (see \code{staircase_steps} parameter). When \code{NULL}
#'  (default), assuming that "\code{adj}" is used as a replacement for
#'  \code{NA}s, \code{adj_init} is calculated as the global alpha divided by the
#'  maximum number of looks (Bonferroni correction), as a rough initial
#'  approximation. However, multiplication is assumed when finding any
#'  multiplication sign (\code{*}) in a given custom \code{adjust} function: in
#'  such a case, \code{adj_init} will be \code{1} by default.
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
#'  any multiplication sign (\code{*}) in a given custom \code{adjust} function.
#'  The former is chosen in any other case.
#'@param alpha_precision During the staircase procedure, at any point when the
#'  simulated global type 1 error rate first matches the given
#'  \code{alpha_global} at least for the number of fractional digits given here
#'  (\code{alpha_precision}; default: \code{5}), the procedure stops and the
#'  results are printed. (Otherwise, the procedures finishes only when all steps
#'  given as \code{staircase_steps} have been used.)
#'@param fut_locals Specifies local futility bounds that may stop the experiment
#'  at the given interim looks if the corresponding p value is above the given
#'  futility bound value. When \code{NULL} (default), sets no futility bounds.
#'  Otherwise, it follows the same logic as \code{alpha_locals} and has the same
#'  input possibilities (number, numeric vector, or named list of numeric
#'  vectors).
#'@param multi_logic_a When multiple p values are evaluated for local alpha
#'  stopping rules, \code{multi_logic_a} specifies the function used for how to
#'  evaluate the multiple significance outcomes (p values being below or above
#'  the given local alphas) as a single \code{TRUE} or \code{FALSE} value that
#'  decides whether or not to stop at a given look. The default, \code{'all'},
#'  specifies that all of the p values must be below the local boundary for
#'  stopping. The other acceptable character input is \code{'any'}, which
#'  specifies that the collection stops when any of the p values pass the
#'  boundary for stopping. Instead of these strings, the actual
#'  \code{\link{all}} and \code{\link{any}} would lead to identical outcomes,
#'  respectively, but the processing would be far slower (since the string
#'  \code{'all'} or \code{'any'} inputs specify a dedicated faster internal
#'  solution). For custom combinations, any custom function can be given, which
#'  will take, as arguments, the p value columns in their given order (either in
#'  the \code{p_values} data frame, or as specified in \code{alpha_locals}), and
#'  should return a single \code{TRUE} or \code{FALSE} value.
#'@param multi_logic_fut Same as \code{multi_logic_a} (again with \code{'all'}
#'  as default), but for futility bounds (for the columns specified in
#'  \code{fut_locals}).
#'@param multi_logic_global Similar as \code{multi_logic_a}, but for the
#'  calculation of the global type 1 error rate (again: in case of multiple p
#'  values being evaluated; otherwise this parameter is not relevant), and with
#'  \code{'any'} by default. This default means that if any of the p values
#'  under evaluation (specified via \code{alpha_locals} or detected
#'  automatically) is significant (p value below the given local alpha at the
#'  stopping of the simulated "experiment" iteration) in case of the H0
#'  scenario, this is calculated as a type 1 error. If \code{'all'} were
#'  specified, only cases with all p evaluated values being significant are
#'  counted as type 1 errors. In either case, the ratio of outcomes with such
#'  type 1 errors (out of all iterations) gives the global type 1 error rate,
#'  which is intended to (approximately) match (via the adjustment procedure)
#'  the value specified as \code{alpha_global}. This global type 1 error is also
#'  what is printed to the console in the end as the "combined" global error
#'  rate. Furthermore, the logic given here is also used for the calculation of
#'  the "combined" global power printed to the console. In this case, the
#'  \code{'any'} logic, for example, would mean that, if any of the p values are
#'  significant at the end of the experiment, this is a positive finding. The
#'  global power is then the ratio of iterations with such positive findings.
#'@param group_by When given as a character element or vector, specifies the
#'  factors by which to group the analysis: the \code{p_values} data will be
#'  divided into parts by these factors and these parts will be analyzed
#'  separately, with power and error information calculated per each part. By
#'  default (\code{NULL}), it identifies factors, if any, given to the
#'  \code{sim} function (via \code{fun_obs}) that created the given
#'  \code{p_values} data.
#'@param alpha_loc_nonstop Optional "non-stopper" alphas via which to evaluate p
#'  values per look, but without stopping the data collection regardless of
#'  statistical significance. Must be a list with names indicating p value
#'  column name pairs, similarly as for the \code{alpha_locals} argument; see
#'  \code{alpha_locals} for details.
#'@param round_to Number of fractional digits (default: \code{5}) to round to,
#'  for the displayed numeric information (such as alphas and power; mainly for
#'  default value for \code{\link[POSSA:print.possa_pow_list]{printing}}).
#'@param iter_limit In some specific cases of unideal/wrong input, the staircase
#'  may get stuck at a given step's loop process. The \code{iter_limit}
#'  parameter specifies the number (by default \code{100}) at which the script
#'  pauses the loop and offers to the user that the procedure be ceased. If the
#'  user chooses to continue, the offer will always be posed again after the
#'  same number of iterations (e.g., by default, after \code{100}, at
#'  \code{200}, then \code{300}, etc.).
#'@param seed Number for \code{\link{set.seed}}; \code{8} by default. Set to
#'  \code{NULL} for random seed.
#'@param hush Logical. If \code{TRUE}, prevents printing any details (or the
#'  progress bar) to console.
#'
#'@return The returns a \code{\link{list}} (with class \code{"possa_pow_list"})
#'  that includes all details of the calculated power, T1ER, and sample
#'  information. This list can be printed legibly (via POSSA's
#'  \code{\link[POSSA:print.possa_pow_list]{print()}} method).
#'
#'@note
#'
#'For the replicability, in case the \code{adjust} function uses any
#'randomization, \code{\link{set.seed}} is executed in the beginning of this
#'function, each time it is called; see the \code{seed} parameter.
#'
#'This function uses, internally, the \code{\link{data.table}} R package.
#'
#'@seealso \code{\link{sim}}
#' @examples
#'
#'# below is a (very) minimal example
#'# for more, see the vignettes via https://github.com/gasparl/possa#usage
#'
#' # create sampling function
#' customSample = function(sampleSize) {
#'     list(
#'         sample1 = rnorm(sampleSize, mean = 0, sd = 10),
#'         sample2_h0 = rnorm(sampleSize, mean = 0, sd = 10),
#'         sample2_h1 = rnorm(sampleSize, mean = 5, sd = 10)
#'     )
#' }
#'
#' # create testing function
#' customTest = function(sample1, sample2_h0, sample2_h1) {
#'  c(
#'    p_h0 = t.test(sample1, sample2_h0, 'less', var.equal = TRUE)$p.value,
#'    p_h1 = t.test(sample1, sample2_h1, 'less', var.equal = TRUE)$p.value
#'  )
#' }
#'
#' # run simulation
#' dfPvals = sim(
#'     fun_obs = customSample,
#'     n_obs = 80,
#'     fun_test = customTest,
#'     n_iter = 1000
#' )
#'
#' # get power info
#' pow(dfPvals)
#'
#'@export
pow = function(p_values,
               alpha_locals = NULL,
               alpha_global = 0.05,
               adjust = TRUE,
               adj_init = NULL,
               staircase_steps = NULL,
               alpha_precision = 5,
               fut_locals = NULL,
               multi_logic_a = 'all',
               multi_logic_fut = 'all',
               multi_logic_global = 'any',
               group_by = NULL,
               alpha_loc_nonstop = NULL,
               round_to = 5,
               iter_limit = 100,
               seed = 8,
               hush = FALSE) {
    all_args = mget(names(formals()), sys.frame(sys.nframe()))
    p_values = data.table::copy(p_values)
    validate_args(
        match.call(),
        list(
            val_arg(p_values, c('df')),
            val_arg(alpha_global, c('num'), 1),
            val_arg(adjust, c('function', 'bool')),
            val_arg(adj_init, c('null', 'num'), 1),
            val_arg(group_by, c('null', 'char')),
            val_arg(multi_logic_a, c('char', 'function'), 1, c('all', 'any')),
            val_arg(multi_logic_fut, c('char', 'function'), 1, c('all', 'any')),
            val_arg(
                multi_logic_global,
                c('char', 'function'),
                1,
                c('all', 'any')
            ),
            val_arg(staircase_steps, c('null', 'num')),
            val_arg(alpha_precision, c('num'), 1),
            val_arg(round_to, c('num'), 1),
            val_arg(iter_limit, c('num'), 1),
            val_arg(seed, c('num'), 1),
            val_arg(hush, c('bool'), 1)
        )
    )
    .look = NULL
    .iter = NULL
    .n_total = NULL
    h0_stoP = NULL
    h1_stoP = NULL
    h0_stoP_sign = NULL
    h1_stoP_sign = NULL
    h0_stoP_fa = NULL
    h1_stoP_fa = NULL
    min_look = NULL
    ._possa_fact_combs = NULL
    . = NULL

    set.seed(seed)
    if (!'possa_sim_df' %in% class(p_values)) {
        warning(
            'The given data frame seems not to have been created by the ',
            '"possa::sim_pvals()" function; it may not fit the "possa::get_pow()" function.',
            immediate. = TRUE
        )
    }
    m_l_reduce = FALSE
    m_l_fut_reduce = FALSE
    if (is.function(multi_logic_a)) {
        if (hush == FALSE) {
            if (isTRUE(all.equal(multi_logic_a, any))) {
                message(
                    'Warning: Operation is hugely faster with ',
                    '"any" string argument for multi_logic_a.'
                )
            } else if (isTRUE(all.equal(multi_logic_a, all))) {
                message(
                    'Warning: Operation is hugely faster with ',
                    '"all" string argument for multi_logic_a.'
                )
            }
        }
    } else {
        if (multi_logic_a == 'all') {
            multi_logic_a = `&`
        } else if (multi_logic_a == 'any') {
            multi_logic_a = `|`
        }
        m_l_reduce = TRUE
    }
    if (is.function(multi_logic_fut)) {
        if (hush == FALSE) {
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
        }
    } else {
        if (multi_logic_fut == 'all') {
            multi_logic_fut = `&`
        } else if (multi_logic_fut == 'any') {
            multi_logic_fut = `|`
        }
        m_l_fut_reduce = TRUE
    }
    m_l_glob_reduce = FALSE
    if (is.function(multi_logic_global)) {
        if (hush == FALSE) {
            if (isTRUE(all.equal(multi_logic_global, any))) {
                message(
                    'Warning: Operation is hugely faster with ',
                    '"any" string argument for multi_logic_global'
                )
            } else if (isTRUE(all.equal(multi_logic_global, all))) {
                message(
                    'Warning: Operation is hugely faster with ',
                    '"all" string argument for multi_logic_global'
                )
            }
        }
    } else {
        if (multi_logic_global == 'all') {
            multi_logic_global = `&`
        } else if (multi_logic_global == 'any') {
            multi_logic_global = `|`
        }
        m_l_glob_reduce = TRUE
    }
    # this is where all outcome dfs will be stored
    out_dfs = list()
    # set data.table
    setDT(p_values)
    setkey(p_values, .look)
    setindex(p_values, .iter)
    looks = unique(p_values$.look)
    mlook = max(p_values$.look)
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
            p_nam = substr(c_nam, 1, nchar(c_nam) - 3)
            if (paste0(p_nam, '_h1') %in% colnames(p_values)) {
                p_names_auto = c(p_names_auto, p_nam)
            }
        }
    }
    a_locals = list()
    p_names = p_names_auto
    # extract (if given) predetermined local alphas and/or specified p value columns
    if (!is.null(alpha_locals)) {
        if (is.atomic(alpha_locals)) {
            # if vector given, assign to each p column
            if (length(alpha_locals) == 1) {
                for (p_nam in p_names) {
                    a_locals[[p_nam]] = rep(alpha_locals, mlook)
                }
            } else if (!length(alpha_locals) == mlook) {
                stop(
                    'Wrong argument for "alpha_locals". (If a numeric vector is given, ',
                    'it must have same length as the maximum number of looks (in this case ',
                    mlook,
                    ').)'
                )
            } else {
                for (p_nam in p_names) {
                    a_locals[[p_nam]] = alpha_locals
                }
            }
            loc_pnames = p_names
        } else {
            # if list, check each vector
            for (a_vec_name in names(alpha_locals)) {
                wrongmes = paste0(
                    'Wrong argument for "alpha_locals". (If a list is given, ',
                    'it must consist of one or more vectors of numbers with the',
                    ' either one value or the same length as the maximum number ',
                    'of looks (in this case ',
                    mlook,
                    ').)'
                )
                if (!(is.atomic(alpha_locals[[a_vec_name]]))) {
                    stop(wrongmes)
                } else if (length(alpha_locals[[a_vec_name]]) == 1) {
                    # if just one value, extend to all looks
                    a_locals[[a_vec_name]] = rep(alpha_locals[[a_vec_name]], mlook)
                } else if (length(alpha_locals[[a_vec_name]]) == mlook) {
                    # if all values, simply assign
                    a_locals[[a_vec_name]] = alpha_locals[[a_vec_name]]
                } else {
                    stop(wrongmes)
                }
            }
            loc_pnames = names(alpha_locals)
        }
        for (pname in loc_pnames) {
            if (!pname %in% p_names) {
                stop(
                    'Wrong argument for "alpha_locals". ',
                    'There is no column name pair "',
                    pname,
                    '_h0"/"',
                    pname,
                    '_h1".'
                )
            }
        }
        p_names = loc_pnames
    }
    if (!length(a_locals) > 0) {
        # if not given, check only last look (with alpha_global)
        # for the rest, local alpha will be 0
        for (p_nam in p_names) {
            a_locals[[p_nam]] = c(rep(0, (mlook - 1)), alpha_global)
        }
        if (isTRUE(adjust)) {
            adjust = FALSE
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
    if (!is.null(alpha_loc_nonstop)) {
        for (vec in alpha_loc_nonstop) {
            if (any(is.na(vec))) {
                stop(
                    'The alpha values given in alpha_loc_nonstop must ',
                    'not contain NA values.'
                )
            } else if (any(vec > 1 | vec < 0)) {
                stop('All alpha values given in alpha_loc_nonstop must be',
                     ' between 0 and 1.')
            } else if (!(is.atomic(vec) &&
                         (length(vec) == mlook ||
                          length(vec) == 1))) {
                stop(
                    'Wrong argument for "alpha_loc_nonstop". (If a list is given, ',
                    'it must consist of one or more vectors of numbers with ',
                    'either one value or the same length as the maximum number ',
                    'of looks (in this case ',
                    mlook,
                    ').)'
                )
            }
            # if only one value is given, extend this value for all looks
            for (ns_vec_name in names(alpha_loc_nonstop)) {
                if (length(alpha_loc_nonstop[[ns_vec_name]]) == 1) {
                    alpha_loc_nonstop[[ns_vec_name]] = rep(alpha_loc_nonstop[[ns_vec_name]], mlook)
                }
            }
        }
        p_extr = names(alpha_loc_nonstop)
        for (pname in p_extr) {
            if (!pname %in% p_names_auto) {
                stop(
                    'Wrong argument for "alpha_loc_nonstop". ',
                    'There is no column name pair "',
                    pname,
                    '0"/"',
                    pname,
                    '1".'
                )
            }
        }
        a_locals = append(a_locals, alpha_loc_nonstop)
    }
    p_names_extr = c(p_names, p_extr)
    fa_locals = list()
    fa_p_names = c()
    # extract (if given) predetermined futility bounds
    if (!is.null(fut_locals)) {
        if (is.atomic(fut_locals)) {
            # if vector given, assign to each p column
            if (length(fut_locals) == 1) {
                for (p_nam in p_names) {
                    fa_locals[[p_nam]] = rep(fut_locals, (mlook - 1))
                }
            } else if (!length(fut_locals) == (mlook - 1)) {
                stop(
                    'Wrong argument for "fut_locals". (If a numeric vector is given, ',
                    'its length must be one less than the maximum number of looks (in this case ',
                    (mlook - 1),
                    ').)'
                )
            } else {
                for (p_nam in p_names) {
                    fa_locals[[p_nam]] = fut_locals
                }
            }
            fa_p_names = p_names
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
                fa_p_names = c(fa_p_names, pname)
            }
            for (p_nam in p_names) {
                if (!p_nam %in% names(fut_locals)) {
                    fut_locals[[p_nam]] = rep(1, (mlook - 1))
                }
            }
            fa_locals = fut_locals
        }
    }
    if (hush == FALSE) {
        prog_bar = TRUE
    } else {
        prog_bar = FALSE
    }
    if (isFALSE(adjust)) {
        # if adjust is FALSE, staircase is omitted
        prog_bar = FALSE
        staircase_steps = NA
    } else if (is.null(staircase_steps)) {
        steps_add = 0.01 * (0.5 ** (seq(0, 11, 1)))
        steps_multi = 0.5 * (0.5 ** (seq(0, 11, 1)))
        if (anyNA(unlist(a_locals))) {
            # default for NA replacement: 11 steps from 0.01, decreasing by halves
            # check: formatC(staircase_steps, digits = 12, format = "f")
            staircase_steps = steps_add
        } else if (isTRUE(adjust)) {
            # given no NAs, assumes multiplication desired
            staircase_steps = steps_multi
            adjust = function(adj, prev, orig) {
                return(orig * adj)
            }
            if (is.null(adj_init)) {
                adj_init = 1 # assumes multiplication
            }
        } else if (is.function(adjust)) {
            staircase_steps = steps_add
            if (any(grepl('*', deparse(body(adjust)), fixed = TRUE))) {
                if (is.null(adj_init)) {
                    adj_init = 1 # assumes multiplication
                }
                # slightly larger steps, again assuming multiplications
                staircase_steps = steps_multi
            }
        } else {
            stop('Wrong argument for "staircase_steps"; please see ?pow.')
        }
    } else if (any(is.na(staircase_steps))) {
        stop('The staircase_steps vector must not contain NA values.')
    }

    if (is.null(adj_init)) {
        # if not given, start adjustment with bonferroni correction
        # (assuming it's going to be used as an alpha level, replacing NAs)
        adj_init = alpha_global / mlook
    }
    if (isTRUE(adjust)) {
        # NA-replacing function
        adjust = function(adj, prev, orig) {
            prev[is.na(orig)] = adj
            return(prev)
        }
    } else if (is.function(adjust)) {
        adjust_args = methods::formalArgs(adjust)
        if (!'adj' %in% adjust_args) {
            stop('The "adjust" function must contain an "adj" parameter.')
        }
        for (a_arg in c('prev', 'orig')) {
            if (!a_arg %in% adjust_args) {
                formals(adjust)[[a_arg]] = NA
            }
        }
        a_example = a_locals[[p_names[1]]]
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
    } else {
        adjust = function(adj, prev, orig) {
            return(orig)
        }
    }

    if (is.null(group_by)) {
        group_by = fac_cols
    } else if (!identical(sort(group_by), sort(fac_cols)) &&
               hush == FALSE) {
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
    p_values[, h0_stoP := TRUE]
    p_values[, h1_stoP := TRUE]
    p_values[, h0_stoP_fa := TRUE]
    p_values[, h1_stoP_fa := TRUE]
    # calculate results separately for each factor combination
    for (possa_fact in possafacts) {
        # when none: possa_fact = NA
        if (is.na(possafacts[1])) {
            pvals_df = p_values
            possa_fact = NULL
        } else {
            # if applicable, take only given factor combination & print its "group" name
            pvals_df = p_values[._possa_fact_combs == possa_fact]
        }
        tot_samples = pvals_df[, stats::median(.n_total), by = '.look']$V1

        ## Sequential design calculation below

        # "trial and error" staircase procedure below, to get the desired adjusted alpha
        locls_temp = a_locals
        a_adj = adj_init
        stair_steps = staircase_steps
        if (!is.null(alpha_loc_nonstop) &
            !is.na(stair_steps[length(stair_steps)])) {
            stair_steps = c(stair_steps, NA)
        }
        a_step = stair_steps[1] # initial a_adj-changing step
        p_h0_sign_names = paste0(p_names, '_h0_sign')
        multi_p = length(p_h0_sign_names) > 1
        p_h1_sign_names = paste0(p_names, '_h1_sign')
        p_h0_fut_names = paste0(fa_p_names, '_h0_fut')
        p_h1_fut_names = paste0(fa_p_names, '_h1_fut')
        if (prog_bar == TRUE) {
            p_bar = utils::txtProgressBar(
                min = 0,
                max = length(staircase_steps),
                initial = 0,
                style = 3
            )
        }
        for (p_nam in p_names) {
            pvals_df[, c(paste0(p_nam, '_h0_sign')) := NA] # create sign column for given p
        }
        for (p_nam in fa_p_names) {
            # if futility bounds are given
            pvals_df[, c(paste0(p_nam, '_h0_fut')) := TRUE] # create fut column for given p
        }
        p_names_temp = p_names
        safe_count = 0
        global_type1 = 1
        while (length(stair_steps) > 0) {
            ### (stepwise info just for testing)
            # cat('\nglobal_type1',
            #     global_type1,
            #     'a_adj:',
            #     a_adj,
            #     'new step:',
            #     a_step,
            #     fill = T)

            # in case of suspiciously many iterations, offer to stop
            safe_count = safe_count + 1
            if (safe_count %% iter_limit == 0) {
                message(
                    paste0(
                        '\nThere have been ',
                        safe_count,
                        ' iterations. The arguments given may be problematic (and/or may never lead to the desired error rate).',
                        '\nNote: The current global type 1 error is ',
                        ro(global_type1, round_to),
                        ' (expected: ',
                        ro(alpha_global, round_to),
                        ').'
                    ),
                    fill = TRUE
                )
                if (readline('Do you want to quit the staircase procedure? (y/n) \n') == 'y') {
                    message('Staircase interrupted.')
                    break
                }
            }

            if (is.na(stair_steps[1])) {
                # as a last step, add non-stopping columns, if any
                p_names_temp = p_names_extr
            } else {
                # adjust local alphas (for stopping alphas only)
                for (p_nam in p_names) {
                    locls_temp[[p_nam]] = adjust(adj = a_adj,
                                                 prev = locls_temp[[p_nam]],
                                                 orig = a_locals[[p_nam]])
                }
            }

            # calculate H0 significances (T/F) & stops (T/F) based on adjusted alphas
            for (p_nam in p_names_temp) {
                for (lk in 1:mlook) {
                    # decide significance at given look for given p
                    pvals_df[.(lk), c(paste0(p_nam, '_h0_sign')) :=
                                 .SD < locls_temp[[p_nam]][lk], .SDcol = paste0(p_nam, '_h0')]
                }
            }
            if (!is.null(fut_locals)) {
                for (p_nam in fa_p_names) {
                    # (mlook - 1 because futility never matters at max look)
                    for (lk in 1:(mlook - 1)) {
                        pvals_df[.(lk), c(paste0(p_nam, '_h0_fut')) :=
                                     .SD > fa_locals[[p_nam]][lk], .SDcol = paste0(p_nam, '_h0')]
                    }
                }
            }

            # check at which look we stop
            # if multiple p values are given as stoppers
            # use multi_logic_a to check where to stop
            # otherwise it simply stops where the given p is sign
            if (multi_p) {
                if (m_l_reduce) {
                    pvals_df[.look != mlook, h0_stoP := Reduce(multi_logic_a, .SD), .SDcols = p_h0_sign_names]
                } else {
                    pvals_df[.look != mlook,  h0_stoP := apply(.SD, 1, multi_logic_a), .SDcols = p_h0_sign_names]
                }
                if (m_l_glob_reduce) {
                    pvals_df[, h0_stoP_sign := Reduce(multi_logic_global, .SD),
                             .SDcols = p_h0_sign_names]
                } else {
                    pvals_df[, h0_stoP_sign := apply(.SD, 1, multi_logic_global),
                             .SDcols = p_h0_sign_names]
                }
            } else {
                pvals_df[.look != mlook, h0_stoP := .SD, .SDcols = p_h0_sign_names]
                pvals_df[, h0_stoP_sign := .SD, .SDcols = p_h0_sign_names]
            }

            if (!is.null(fut_locals)) {
                # analogue with futility bounds
                if (multi_p) {
                    if (m_l_fut_reduce) {
                        pvals_df[.look != mlook, h0_stoP_fa := Reduce(multi_logic_fut, .SD), .SDcols = p_h0_fut_names]
                    } else {
                        pvals_df[.look != mlook, h0_stoP_fa := apply(.SD, 1, multi_logic_fut), .SDcols = p_h0_fut_names]
                    }
                } else {
                    pvals_df[.look != mlook, h0_stoP_fa := .SD, .SDcols = p_h0_fut_names]
                }
                pvals_df[.look != mlook, h0_stoP := h0_stoP |
                             h0_stoP_fa]
            }

            # now get all outcomes at stopping point
            pvals_stp = pvals_df[h0_stoP == TRUE, .SD,
                                 .SDcols = c('h0_stoP_sign', '.look', '.iter')]
            # the global type 1 error
            global_type1 = mean(pvals_stp[, min_look := min(.look), by = .iter][.look == min_look, h0_stoP_sign])

            #mean(pvals_stp[, min_look := min(.look), by = .iter][.look == min_look, p_test2_h0_sign])
            #mean(pvals_stp[, min_look := min(.look), by = .iter][.look == min_look, p_test1_h0_sign])

            if (is.na(stair_steps[1])) {
                # NA indicates no stairs; nothing left to be done
                break
            } else {
                # break if global alpha is correct
                # otherwise continue the staircase
                if (round(global_type1, alpha_precision) == round(alpha_global, alpha_precision)) {
                    if (is.null(alpha_loc_nonstop)) {
                        break
                    } else {
                        stair_steps = NA
                        a_step = 0
                    }
                } else if ((global_type1 < alpha_global &&
                            a_step < 0) ||
                           (global_type1 > alpha_global &&
                            a_step > 0)) {
                    # change staircase direction (and also decrease step) if needed
                    a_step = -stair_steps[1] * sign(a_step)
                    stair_steps = stair_steps[-1]
                    if (hush == FALSE) {
                        utils::setTxtProgressBar(p_bar,
                                                 utils::getTxtProgressBar(p_bar) + 1)
                    }
                }
                # adjust a_adj itself by the step
                # this a_adj will be used in the adjust function
                a_adj = a_adj + a_step
            }
        }
        if (prog_bar == TRUE) {
            utils::setTxtProgressBar(p_bar, length(staircase_steps))
            close(p_bar)
        }
        # assign final local alphas
        a_locals_fin = locls_temp

        for (a_name in names(a_locals_fin)) {
            a_locals_fin[[a_name]][a_locals_fin[[a_name]] < 0] = 0
        }

        # calculate H1 significances (T/F) & stops (T/F) based on final alphas
        # (analogue of the H0 significances above)
        for (p_nam in p_names_extr) {
            pvals_df[, c(paste0(p_nam, '_h1_sign')) := NA]
            for (lk in 1:mlook) {
                pvals_df[.(lk), c(paste0(p_nam, '_h1_sign')) :=
                             .SD < locls_temp[[p_nam]][lk], .SDcol = paste0(p_nam, '_h1')]
            }
        }
        if (!is.null(fut_locals)) {
            for (p_nam in fa_p_names) {
                pvals_df[, c(paste0(p_nam, '_h1_fut')) := TRUE]
                # (mlook - 1 because futility never matters at max look)
                for (lk in 1:(mlook - 1)) {
                    pvals_df[.(lk), c(paste0(p_nam, '_h1_fut')) :=
                                 .SD > fa_locals[[p_nam]][lk], .SDcol = paste0(p_nam, '_h1')]
                }
            }
        }

        # now check the global power
        # if multiple p columns, check at which look we stop
        if (multi_p) {
            if (m_l_reduce) {
                pvals_df[.look != mlook, h1_stoP := Reduce(multi_logic_a, .SD), .SDcols = p_h1_sign_names]
            } else {
                pvals_df[.look != mlook, h1_stoP := apply(.SD, 1, multi_logic_a), .SDcols = p_h1_sign_names]
            }
            if (m_l_glob_reduce) {
                pvals_df[, h1_stoP_sign := Reduce(multi_logic_global, .SD),
                         .SDcols = p_h1_sign_names]
            } else {
                pvals_df[, h1_stoP_sign := apply(.SD, 1, multi_logic_global),
                         .SDcols = p_h1_sign_names]
            }
        } else {
            pvals_df[.look != mlook, h1_stoP := .SD, .SDcols = p_h1_sign_names]
            pvals_df[, h1_stoP_sign := .SD, .SDcols = p_h1_sign_names]
        }

        if (!is.null(fut_locals)) {
            if (multi_p) {
                if (m_l_fut_reduce) {
                    pvals_df[.look != mlook, h1_stoP_fa := Reduce(multi_logic_fut, .SD), .SDcols = p_h1_fut_names]
                } else {
                    pvals_df[.look != mlook,  h1_stoP_fa := apply(.SD, 1, multi_logic_fut), .SDcols = p_h1_fut_names]
                }
            } else {
                pvals_df[.look != mlook, h1_stoP_fa := .SD, .SDcols = p_h1_fut_names]
            }
            pvals_df[.look != mlook, h1_stoP := h1_stoP |
                         h1_stoP_fa]
        }
        if (multi_p) {
            # if multi_p, get global power at stopping point
            pvals_stp = pvals_df[h1_stoP == TRUE, .SD,
                                 .SDcols = c('h1_stoP_sign', '.look', '.iter')]
            # the global type 1 error (not necessary to calculate; will be in totals)
            # global_power = mean(pvals_stp[, min_look := min(.look), by = .iter][.look == min_look, h1_stoP_sign])
        }

        # calculate sample size information per look
        ps_sub0 = data.table::copy(pvals_df)
        ps_sub1 = data.table::copy(pvals_df)
        iters_tot = length(unique(pvals_df$.iter))
        stops = list() # collect info per each stop
        previous_h0 = iters_tot # start with max for both
        previous_h1 = iters_tot
        for (lk in looks) {
            # get iterations stopped at given look
            iters_out0 = ps_sub0[.look == lk &
                                     h0_stoP == TRUE]
            # remove stopped iterations
            ps_sub0 = ps_sub0[!.iter %in% iters_out0$.iter, ]
            # (same for H1)
            iters_out1 = ps_sub1[.look == lk &
                                     h1_stoP == TRUE,]
            ps_sub1 = ps_sub1[!.iter %in% iters_out1$.iter,]
            outs = c()
            # get combined significants
            outs['iters_sign_h0'] = nrow(iters_out0[h0_stoP_sign == TRUE])
            outs['iters_sign_h1'] = nrow(iters_out1[h1_stoP_sign == TRUE])
            if (!is.null(fut_locals)) {
                outs['iters_fut_h0'] = nrow(iters_out0[h0_stoP_fa == TRUE])
                outs['iters_fut_h1'] = nrow(iters_out1[h1_stoP_fa == TRUE])
            }
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
                length(unique(ps_sub0$.iter))
            outs['iters_remain_h1'] =
                length(unique(ps_sub1$.iter))
            if (lk == mlook) {
                # at last look, all is stopped that previously remained
                outs['iters_stopped_h0'] = previous_h0
                outs['iters_stopped_h1'] = previous_h1
            } else {
                # calculate stops as previous minus current remaining
                outs['iters_stopped_h0'] = previous_h0 - outs['iters_remain_h0']
                outs['iters_stopped_h1'] = previous_h1 - outs['iters_remain_h1']
            }
            stops[[length(stops) + 1]] = c(
                look = lk,
                n = tot_samples[lk],
                n_rate = tot_samples[lk] / tot_samples[mlook],
                outs
            )
            # assign current remaining as the next "previous remaining"
            previous_h0 = outs['iters_remain_h0']
            previous_h1 = outs['iters_remain_h1']
        }
        df_stops = as.data.frame(do.call(rbind, stops))
        # derive end info (type 1, power, etc) from the iteration ratios
        for (p_nam in p_names_extr) {
            # write out local alphas
            df_stops[paste0('alpha_local_', p_nam)] = a_locals_fin[[p_nam]]
            if (!is.null(fut_locals)) {
                # write out local futility bounds (if any)
                df_stops[paste0('futil_local_', p_nam)] = c(fa_locals[[p_nam]], NA)
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

        # count ratios of combined stopping significances per each look
        df_stops$ratio_combined_sign_h0 = df_stops$iters_sign_h0 / iters_tot
        df_stops$ratio_combined_sign_h1 = df_stops$iters_sign_h1 / iters_tot

        # count ratios of combined futility stops (if any) per each look
        if (!is.null(fut_locals)) {
            df_stops$ratio_combined_fut_h0 = df_stops$iters_fut_h0 / iters_tot
            df_stops$ratio_combined_fut_h1 = df_stops$iters_fut_h1 / iters_tot
        }

        # count cumulative ratios of remainings
        df_stops$ratio_remain_h0 = df_stops$iters_remain_h0 / iters_tot
        df_stops$ratio_remain_h1 = df_stops$iters_remain_h1 / iters_tot
        # the average sample proportion at the given look (mainly just to calculate average-total)
        df_stops$n_avg_prop_0 = df_stops$ratio_stopped_h0 * df_stops$n
        df_stops$n_avg_prop_1 = df_stops$ratio_stopped_h1 * df_stops$n
        # calculate sums per each column (with different meaning in each case)
        df_stops = rbind(df_stops, colSums(df_stops))
        df_stops$look[nrow(df_stops)] = 'totals'

        for (p_nam in p_names_extr) {
            if (p_nam %in% p_names) {
                class(df_stops[[paste0('alpha_local_',
                                       p_nam)]]) = c(class(df_stops[[paste0('alpha_local_', p_nam)]]),
                                                     'possa_p')
            } else {
                class(df_stops[[paste0('alpha_local_',
                                       p_nam)]]) = c(class(df_stops[[paste0('alpha_local_', p_nam)]]),
                                                     'possa_p_nonstopper')
            }
            if (!is.null(df_stops[[paste0('futil_local_', p_nam)]])) {
                class(df_stops[[paste0('futil_local_',
                                       p_nam)]]) = c(class(df_stops[[paste0('alpha_local_', p_nam)]]),
                                                     'possa_futility')
            }
        }
        class(df_stops) = c("possa_pow_df", "data.frame")
        if (is.null(possa_fact)) {
            out_dfs$pow = df_stops
        } else {
            out_dfs[[paste('pow', gsub('; ', '_', possa_fact), sep = '_')]] = df_stops
        }
    }
    out_dfs$arguments = all_args
    class(out_dfs) = c("possa_pow_list", class(out_dfs))
    return(out_dfs)
}
