#'@title Simulation procedure
#'
#'@description This function performs the simulation procedure in order to get
#'  the p values that will eventually serve for power calculations (via
#'  \code{\link{pow}}). The observation values ("sample") to be tested are
#'  simulated via the given \code{fun_obs} function, and the significance
#'  testing is performed via the given \code{fun_test} function. The numbers of
#'  observations per look (for a sequential design) are specified in
#'  \code{n_obs}.
#'
#'@param n_obs A numeric vector or a named list of numeric vectors. Specifies
#'  the numbers of observations (i.e., samples sizes) that are to be generated
#'  by \code{fun_obs} and then tested in \code{fun_test}. If a single vector is
#'  given, this will be used for all observation number arguments in the
#'  \code{fun_obs} and for the sample size adjustments for the arguments in the
#'  \code{fun_test} functions. Otherwise, if a named list of numeric vectors is
#'  given, the names must correspond exactly to the argument names in
#'  \code{fun_obs} and \code{fun_test}, so that the respective numeric vectors
#'  are used for each given sample variable. For convenience, in case of a
#'  "\code{_h}" suffix, the variable will be divided into names with
#'  "\code{_h0}" and "\code{_h1}" suffixes for \code{fun_test} (but not for
#'  \code{fun_obs}); see Details.
#'@param fun_obs A \code{\link{function}} that creates the observations (i.e.,
#'  the "sample"; all values for the dependent variable(s)). The respective
#'  maximum observation number(s), given in \code{n_obs}, will be passed to the
#'  \code{fun_obs}. For this, the returned value must be a named list, where the
#'  names correspond exactly to the arguments in \code{fun_test}. In case of
#'  sequential testing, the observations returned by \code{fun_obs} will be
#'  reduced to the specified (smaller) number(s) of observations for each given
#'  interim "look" (as a simulation for what would happen if collection was
#'  stopped at that given look), to be used in \code{fun_test}. Optionally, the
#'  \code{fun_obs} can be passed additional arguments (via a
#'  \code{\link{list}}); see Details.
#'@param fun_test The function for significance testing. The list of samples
#'  returned by \code{fun_obs} (with observation numbers specified in
#'  \code{n_obs}) will be passed into this \code{fun_test} function as
#'  arguments, to be used in the given statistical significance tests in this
#'  function. To correctly calculate the sample sizes in
#'  \code{\link[POSSA:pow]{POSSA::pow}}, the argument names for the sample that
#'  varies depending on whether the null (H0) and alternative (H1) hypothesis is
#'  true should be indicated with "\code{_h0}" and "\code{_h1}" suffixes,
#'  respectively, with a common root (so, e.g., "\code{var_x_h0}" and
#'  "\code{var_x_h1}"). Then, in the resulting \code{\link{data.frame}}, their
#'  sample size (which must always be identical) will be automatically merged
#'  into a single column with a trimmed "\code{_h}" suffix (e.g.,
#'  "\code{var_x_h}"). (Otherwise, the sample sizes of both H0 and H1 would be
#'  calculated toward the total expected sample in either case, which is of
#'  course incorrect. There are internal checks to prevent this, but the
#'  intended total sample size can also be double-checked in the returned
#'  \code{\link{data.frame}}'s \code{.n_total} column.) Within-subject
#'  observations, i.e., multiple observations per group, should be specified
#'  with "\code{GRP}" prefix for a single group (e.g., simply "\code{GRP}", or
#'  "\code{GRP_mytest}") and, for multiple groups, "\code{grp_}" prefix with a
#'  following group name (e.g., "\code{grp_1}" or "\code{grp_alpha}"); the
#'  numbers of multiple observations in each group can then be specified in
#'  \code{fun_obs} via their group name (since the respective numbers of
#'  observations should always be the same anyway); see Examples. To be
#'  recognized by the \code{\link[POSSA:pow]{POSSA::pow}} function, the
#'  \code{fun_test} must return a named vector including a pair (or pairs) of p
#'  values for H0 and H1 outcomes, where each p value's name must be specified
#'  with a "\code{p_}" prefix and a "\code{_h0}" suffix for H0 outcome or a
#'  "\code{_h1}" suffix for H1 outcome (e.g., \code{p_h0}, \code{p_h1};
#'  \code{p_ttest_h0}, \code{p_ttest_h1}). The simulated outcomes (per
#'  iteration) for each of these p values will be separately stored in a
#'  dedicated column of the \code{\link{data.frame}} returned by the \code{sim}
#'  function. Optionally, the \code{fun_test} can return other miscellaneous
#'  outcomes too, such as effect sizes or confidence interval limits; these will
#'  then be stored in dedicated columns in the resulting
#'  \code{\link{data.frame}}.
#'@param n_iter Number of iterations (default: 45000).
#'@param adjust_n Adjust total number of observations via simple multiplication.
#'  Might be useful in some specific cases, e.g. if for some reason multiple p
#'  values are derived from the same sample without specifying grouping
#'  (\code{GRP} or \code{grp_} in \code{fun_test}), which would then lead to
#'  incorrect (too many, multiplied) totals; for example, in case of four
#'  observations obtained from the same sample, the value \code{1/4} could be
#'  given. (The default value is \code{1}.)
#'@param seed Number for \code{\link{set.seed}}; \code{8} by default. Set to
#'  \code{NULL} for random seed.
#'@param pair Logical or \code{NULL}. By default \code{NULL}, the algorithm
#'  assumes paired samples included among the observations in case of any
#'  grouping via the \code{fun_test} parameters ("\code{GRP}"/"\code{grp}"), and
#'  no paired samples otherwise. In case of paired samples included, within each
#'  look, the same vector indexes to remove elements from the given
#'  observations. In general, this should not substantially affect the outcomes
#'  of independent samples (assuming that their order is truly independent), but
#'  this depends on how the random samples are generated in the \code{fun_obs}
#'  function. To be safe and avoid any potential bias, it is best to avoid this
#'  paired sampling mechanism when no paired samples are included. To override
#'  the default, set to \code{TRUE} for paired samples scenario (paired
#'  sampling), or to \code{FALSE} for no paired samples scenario (random
#'  subsampling of each sample). (Might be useful for testing or some very
#'  specific procedures, e.g. where grouping is not indicated despite paired
#'  samples.)
#'@param ignore_suffix Set to \code{NULL} to give warnings instead of errors for
#'  internally detected consistency problems with the \code{_h0}/\code{_h1}
#'  suffixes in the \code{fun_test} function arguments. Set to \code{TRUE} to
#'  completely ignore these (neither error nor warning). (Might be useful for
#'  testing or some very specific procedures.)
#'@param prog_bar Logical, \code{FALSE} by default. If \code{TRUE}, shows
#'  progress bar.
#'@param hush Logical, \code{FALSE} by default. If \code{TRUE}, prevents
#'  printing any details to console.
#'
#'@details
#'
#'To specify a variable that differs depending on whether the null hypothesis
#'("H0") or the alternative hypothesis ("H1") is true, a pair of samples are
#'needed for \code{fun_test}, for which the argument names should have an
#'identical root and "\code{_h0}" and "\code{_h1}" endings, such as
#'"\code{var_x_h0}" (for sample in case of H0) and "\code{var_x_h1}" (for sample
#'in case of H1). Then, since the observation number for this pair will always
#'be the same, as a convenience, parameters with "\code{_h0}" and "\code{_h1}"
#'endings specifically can be specified together in \code{n_obs} with the last
#'"0"/"1" character dropped, hence ending with "\code{_h}". So, for example,
#'"\code{var_x_h = c(30, 60, 90)}" will be automatically adjusted to specify the
#'observation numbers for both "\code{var_x_h0}" and "\code{var_x_h1}". In that
#'case, \code{fun_obs} must have a single argument "\code{var_x_h}", while
#'\code{fun_test} must have both full names as arguments ("\code{var_x_h0}" and
#'"\code{var_x_h1}").
#'
#'Optionally, \code{fun_obs} can be provided in \code{\link{list}} format for
#'the convenience of exploring varying factors (e.g., different effect sizes,
#'correlations) at once, without writing a dedicated \code{fun_obs} function for
#'each combination, and each time separately running the simulation and the
#'power calculation. In this case, the first element of the list must be the
#'actual \code{\link{function}}, which contains certain parameters for
#'specifying varying factors, while the rest of the elements should contain the
#'various argument values for these parameters of the function as named elements
#'of the list (e.g., \code{list(my_function, factor1=c(1, 2, 3), factor2=c(0,
#'5))}), with the name corresponding to the parameter name in the function, and
#'the varying values (numbers or strings). When so specified, a separate
#'simulation procedure will be run for each combination of the given factors
#'(or, if only one factor is given, for each element of that factor). The
#'\code{\link[POSSA:pow]{POSSA::pow}} function will be able to automatically
#'detect (by default) the factors generated this way in the present
#'\code{\link[POSSA:sim]{POSSA::sim}} function, in order to calculate power
#'separately for each factor combination.
#'
#'@return Returns a \code{\link{data.frame}} (with class \code{"possa_sim_df"})
#'  that includes the columns \code{.iter} (the iterations of the simulation
#'  procedure numbered from \code{1} to \code{n_iter}), \code{.look} (the
#'  interim "looks" numbered from \code{1} to the maximum number of looks,
#'  including the final one), and the information returned by the
#'  \code{fun_test} function for H0 and H1 outcomes (mainly p values; but also
#'  other, optional information, if any) and the corresponding observation
#'  numbers, as well as the total observation number per each look under a
#'  dedicated \code{.n_total} column. When this data frame is printed to the
#'  console (via POSSA's \code{\link[POSSA:print.possa_pow_list]{print()}}
#'  method), the head (first few lines) of the data is shown, as well as, in
#'  case of any varying factors included, summary information per factor
#'  combination.
#'
#'@note
#'
#'For the replicability (despite the randomization), \code{\link{set.seed}} is
#'executed in the beginning of this function, each time it is called; see the
#'\code{seed} parameter.
#'
#'@seealso \code{\link{pow}}
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
#'
#' @export
sim = function(fun_obs,
               n_obs,
               fun_test,
               n_iter = 45000,
               adjust_n = 1,
               seed = 8,
               pair = NULL,
               ignore_suffix = FALSE,
               prog_bar = FALSE,
               hush = FALSE) {
    validate_args(
        match.call(),
        list(
            val_arg(fun_obs, c('function', 'list')),
            val_arg(fun_test, c('function')),
            val_arg(n_iter, c('num'), 1),
            val_arg(adjust_n, c('num'), 1),
            val_arg(seed, c('null', 'num'), 1),
            val_arg(pair, c('bool', 'null'), 1),
            val_arg(ignore_suffix, c('null', 'bool'), 1),
            val_arg(prog_bar, c('bool'), 1),
            val_arg(hush, c('bool'), 1)
        )
    )
    set.seed(seed)
    # sanity checks for given values
    if (!is.function(fun_obs)) {
        f_obs_args = fun_obs
        if (!is.function(f_obs_args[[1]])) {
            stop('When "fun_obs" is given as list,',
                 ' the first argument must be a function.')
        }
        fun_obs = f_obs_args[[1]]
        f_obs_args[[1]] = NULL
        # get all fun_obs argument combinations (different sample types)
        df_combs = sapply(expand.grid(f_obs_args), as.vector)
        facts_list = list()
        for (rownum in 1:nrow(df_combs)) {
            facts_list[[rownum]] = as.list(df_combs[rownum,])
        }
    } else {
        # set to have no combinations; single sample test (hence 1 cycle below)
        facts_list = NA
        f_obs_args = c()
    }
    f_test_Arg_names = methods::formalArgs(fun_test)
    f_obs_Arg_names = methods::formalArgs(fun_obs)
    f_obs_Arg_names = f_obs_Arg_names[!f_obs_Arg_names %in% names(f_obs_args)]
    # get groups based on f_test and f_obs arguments
    grp_samples = list()
    grp_obs = list()
    for (arg_name in f_test_Arg_names) {
        # if any "grp" name is given, assign same numbers for each given group
        if (startsWith(arg_name, 'grp_') |
            startsWith(arg_name, 'GRP')) {
            if (startsWith(arg_name, 'grp_')) {
                grp_name = paste0('grp_', strsplit(arg_name, '_')[[1]][2])
            } else {
                grp_name = 'GRP'
            }
            # if group not yet created, create now
            if (!grp_name %in% names(grp_samples)) {
                # all corresponding group sample names guessed from all fun_test arguments
                grp_Arg_Test = f_test_Arg_names[startsWith(f_test_Arg_names, grp_name)]
                if (length(grp_Arg_Test) > 0) {
                    grp_samples[[grp_name]] = grp_Arg_Test
                    if (hush == FALSE) {
                        message(
                            'Note: Observation numbers groupped as "',
                            grp_name,
                            '" for ',
                            paste(grp_Arg_Test, collapse = ', '),
                            '.'
                        )
                    }
                    if (is.null(pair)) {
                        pair = TRUE
                    }
                }
            }
        }
    }
    if (is.null(pair)) {
        pair = FALSE
    }
    n_obs_max = list()
    if (is.atomic(n_obs)) {
        n_look = length(n_obs)
        # if just a vector given, all samples have this as samples sizes
        n_obs_orig = n_obs
        n_obs = list()
        # required sample names guessed from fun_test arguments
        for (n_name in f_test_Arg_names) {
            n_obs[[n_name]] = n_obs_orig # assign vector to each sample type
        }
        # required sample size names guessed from fun_obs arguments
        for (n_name in f_obs_Arg_names) {
            if (!n_name %in% names(f_obs_args)) {
                n_obs_max[[n_name]] = n_obs_orig[n_look] # assign vector to each sample type
            }
        }
    } else {
        n_look = length(n_obs[[1]])
        for (n_name in names(n_obs)) {
            n_obs_max[[n_name]] = n_obs[[n_name]][n_look]
            if (n_name %in% names(grp_samples)) {
                # if grouped, assign each subgroup for n_obs_max
                for (subgrp in grp_samples[[n_name]]) {
                    n_obs[[subgrp]] = n_obs[[n_name]]
                }
            }
            if (endsWith(n_name, '_h')) {
                n_obs[[paste0(n_name, '0')]] = n_obs[[n_name]]
                n_obs[[paste0(n_name, '1')]] = n_obs[[n_name]]
            }
        }
        if (!identical(sort(names(n_obs_max)), sort(f_obs_Arg_names))) {
            stop(
                'The names provided via "n_obs" (',
                paste(names(n_obs_max), collapse = ', '),
                ') do not match the arguments in "fun_obs" (',
                paste(f_obs_Arg_names, collapse = ', '),
                ').'
            )
        }
        # remove unnecessary duplicates in n_obs, depending on what fun_test needs
        n_obs = n_obs[names(n_obs) %in% f_test_Arg_names]
        if (!identical(sort(names(n_obs)), sort(f_test_Arg_names))) {
            stop(
                'The names provided via "n_obs" (',
                paste(names(n_obs), collapse = ', '),
                ') do not match the arguments in "fun_test" (',
                paste(f_test_Arg_names, collapse = ', '),
                ').'
            )
        }
    }
    if (is.null(ignore_suffix)) {
        feedf = warning
    } else if (isFALSE(ignore_suffix)) {
        feedf = stop
    } else {
        feedf = function(...) {
            # nothing
        }
    }
    obs_root_names = suffixes(f_test_Arg_names,
                              feedf = feedf)
    if (!length(obs_root_names) > 0) {
        feedf(
            'No "_h0"/"_h1" suffix pair found among ',
            'the "fun_test" arguments. Please see ?sim."'
        )
    }
    if (!is.na(facts_list[1])) {
        fun_test_out = do.call(fun_test, do.call(fun_obs, c(n_obs_max, facts_list[[1]])))
    } else {
        fun_test_out = do.call(fun_test, do.call(fun_obs, c(n_obs_max)))
    }
    p_root_names = suffixes(names(fun_test_out)[startsWith(names(fun_test_out), 'p_')], feedf = feedf)
    if (!(is.vector(fun_test_out) &
          !is.null(names(fun_test_out)) &
          !any(is.na(names(fun_test_out))))) {
        feedf('The "fun_test" function must return a named vector.',
              ' See ?sim for details."')
    }
    if (!length(p_root_names) > 0) {
        feedf(
            'No "_h0"/"_h1" suffix pair found among ',
            'the "fun_test" results. Please see ?sim."'
        )
    }

    if (!all(sapply(n_obs, length) == n_look)) {
        stop('The lengths of the "n_obs" values are unequal.')
    }
    obs_per_it = do.call(Map, c(f = list, n_obs)) # transpose input for iteration below
    obs_names = names(obs_per_it[[1]]) # the names of all sample columns
    list_vals = list() # all end results will be collected in this list
    
    if (hush == TRUE) {
        prog_bar = FALSE
    }
    
    # set progress bar
    if (prog_bar == TRUE) {
        pb = utils::txtProgressBar(
            min = 0,
            max = n_iter * length(facts_list),
            initial = 0,
            style = 3
        )
    }
    pb_count = 0
    if (n_look > 1) {
        looks = (n_look - 1):1
    } else {
        looks = NULL
    }

    # pre-compile main functions
    # (doesn't seem to actually make any difference, but can't hurt)
    fun_obs = compiler::cmpfun(fun_obs)
    fun_test = compiler::cmpfun(fun_test)

    # start power calculation per each factor combination
    for (facts in facts_list) {
        if (is.na(facts[1])) {
            facts = NULL
        }
        for (i in 1:n_iter) {
            pb_count = pb_count + 1
            if (prog_bar == TRUE) {
                utils::setTxtProgressBar(pb, pb_count)
            }
            # call the full (max) samples
            samples = do.call(fun_obs, c(n_obs_max, facts))
            # test the full samples
            list_vals[[length(list_vals) + 1]] =
                c(
                    .iter = i,
                    .look = n_look,
                    unlist(obs_per_it[[n_look]]),
                    do.call(fun_test, samples)
                )
            # now, subsample each full sample, per each look
            for (lk in looks) {
                if (pair == TRUE) {
                    seed_r = .GlobalEnv$.Random.seed # .Random.seed
                }
                for (samp_n in obs_names) {
                    if (pair == TRUE) {
                        .GlobalEnv$.Random.seed = seed_r # assign(".Random.seed", seed_r, envir = parent.frame())
                    }
                    # take subsample for the given look with the given size
                    samples[[samp_n]] = sample(samples[[samp_n]],
                                               obs_per_it[[lk]][[samp_n]])
                }
                # repeat the test with reduced samples
                list_vals[[length(list_vals) + 1]] =
                    c(
                        .iter = i,
                        .look = lk,
                        unlist(obs_per_it[[lk]]),
                        do.call(fun_test, samples)
                    )
            }
        }
        df_fact = as.data.frame(do.call(rbind, list_vals))
        if (is.null(facts)) {
            df_pvals = df_fact
        } else if (exists("df_pvals")) {
            df_pvals = rbind(df_pvals,
                             data.frame(df_fact[1:2], facts, df_fact[c(-1,-2)]))
        } else {
            df_pvals = data.frame(df_fact[1:2], facts, df_fact[c(-1,-2)])
        }
        list_vals = list()
    }
    if (prog_bar == TRUE) {
        close(pb)
    }

    # merge groups into single column (with given group name)
    for (grp_nam in names(grp_samples)) {
        obs_colnames = grp_samples[[grp_nam]] # get group columns
        # check if all group columns have the same number of observations
        if (!length(unique(as.list(df_pvals[obs_colnames]))) == 1) {
            warning(
                'Observation numbers do not match for the columns ',
                paste(obs_colnames, collapse = ', '),
                ' (group: "',
                grp_nam,
                '")! Merging columns omitted;',
                ' total sample size may be incorrect.'
            )
        } else {
            # renaming the first group member to represent all (since all are identical)
            names(df_pvals)[names(df_pvals) == obs_colnames[1]] = grp_nam
            # remove other group columns from dataframe
            df_pvals = df_pvals[, !(names(df_pvals) %in% obs_colnames)]
            # remove all group columns from obs names
            obs_names = obs_names[!obs_names %in% obs_colnames]
            # add group name to obs names, to represent all group columns
            obs_names = c(obs_names, grp_nam)
        }
    }
    # merge _h0/_h1 into single _h column
    for (obs_nam in obs_root_names) {
        obs_pair = c(paste0(obs_nam, '0'), paste0(obs_nam, '1'))
        if (obs_pair[1] %in% names(df_pvals)) {
            if (all(df_pvals[[obs_pair[1]]] != df_pvals[[obs_pair[2]]])) {
                warning(
                    'Observation numbers do not match for the ',
                    obs_nam,
                    ' h0/h1 column pair! Merging columns omitted;',
                    ' total sample size may be incorrect.'
                )
            } else {
                # renaming the first to represent both (since they are identical)
                names(df_pvals)[names(df_pvals) == obs_pair[1]] = obs_nam
                # removing the other
                df_pvals[[obs_pair[2]]] = NULL
                # removing the pair from obs names
                obs_names = obs_names[!obs_names %in% obs_pair]
                # adding pair root name to represent both
                obs_names = c(obs_names, obs_nam)
            }
        }
    }
    # calculate .n_total (total observation number per each "look")
    if (length(obs_names) > 1) {
        n_tots = Reduce('+', df_pvals[, obs_names])
    } else {
        n_tots =  df_pvals[, obs_names]
    }
    df_pvals = data.frame(df_pvals[, 1:2],
                          .n_total = adjust_n * n_tots,
                          df_pvals[,-1:-2])
    # order per iter and look
    df_pvals = df_pvals[order(df_pvals$.iter, df_pvals$.look),]
    if (length(f_obs_args) > 0) {
        df_pvals = df_pvals[do.call('order', df_pvals[names(f_obs_args)]), ]
    }

    # add POSSA class names, to be recognized in POSSA::pow
    for (c_nam in obs_names) {
        # observation number (sample size) columns
        class(df_pvals[[c_nam]]) = c(class(df_pvals[[c_nam]]), 'possa_n')
    }
    for (fc_nam in names(f_obs_args)) {
        # columns specifying varying factors
        class(df_pvals[[fc_nam]]) = c(class(df_pvals[[fc_nam]]), 'possa_fac')
    }
    # POSSA class for the whole data frame
    class(df_pvals) = c('possa_sim_df', class(df_pvals))
    return(df_pvals)
}