#'@title Simulation procedure
#'
#'@description This function performs the simulation procedure in order to get
#'  the p values that will eventually serve for power calculations (via
#'  \code{\link{pow}}). The observation values ("sample") to be tested are
#'  simulated via the given \code{fun_obs} function, and the significance
#'  testing is performed via the given \code{fun_test} function. The numbers of
#'  observations per look (for a sequential design) are specified in
#'  \code{n_obs}.
#'@param fun_obs A \code{\link{function}} that creates the observations (i.e.,
#'  the "sample"; all values for the dependent variable(s)). The arguments that
#'  specify the observation numbers in this function have to be provided via
#'  \code{n_obs}: the respective maximum observation number, given in
#'  \code{n_obs}, will be passed to the \code{fun_obs}, and, in case of
#'  sequential testing, the resulting observations will be reduced to the
#'  specified (smaller) number(s) of observations for each given interim "look"
#'  (as a simulation for what would happen if collection was stopped at that
#'  given look). Importantly, a special "\code{_h}" suffix (e.g.,
#'  "\code{var1_h}") is required to specify the variable(s) that differ(s)
#'  depending on whether or not the null hypothesis ("H0") or the alternative
#'  hypothesis ("H1") is true. This will be internally transformed into
#'  "\code{_h0}" and "\code{_h1}" endings (e.g., \code{var1_h0}, \code{var1_h1})
#'  for the arguments of the \code{fun_test} function.  Optionally, the
#'  \code{fun_obs} can be passed additional arguments (via a
#'  \code{\link{list}}); see Details.
#'@param n_obs A numeric vector or a (named) list of numeric vectors. Specifies
#'  the numbers of observations (i.e., samples sizes) that are to be generated
#'  from \code{fun_obs}, to be then tested in \code{fun_test}. If a single
#'  vector is given, this will be used for all (observation number) arguments in
#'  the \code{fun_obs} function. Otherwise, if a named list of numeric vectors
#'  is given, the names must correspond exactly to the argument names in
#'  \code{fun_obs}, so that the respective numeric vectors are used for each
#'  given variable.
#'@param fun_test The function for significance testing. The observation values
#'  created in \code{fun_obs} (with observation numbers specified in
#'  \code{n_obs}) will be passed into this \code{fun_test} function as arguments
#'  (with \code{_h} suffixes extended to \code{_h0} and \code{_h1}), to be used
#'  in the given statistical significance tests in this function. The test
#'  results must include a pair (or pairs) of p values for H0 and H1 outcomes,
#'  where each p value must be specified with a "\code{p_}" prefix and a
#'  "\code{_h0}" suffix for H0 outcome or a "\code{_h1}" suffix for H1 outcome
#'  (e.g., \code{p_h0}, \code{p_h1}; \code{p_ttest_h0}, \code{p_ttest_h1}). The
#'  simulated outcomes (per iteration) for each of these p values will be
#'  separately stored in a dedicated column of the \code{\link{data.frame}}
#'  returned by the \code{sim} function. Optionally, the \code{fun_test} can
#'  return other miscellaneous outcomes too, such as effect sizes or confidence
#'  interval limits; these will then be stored in dedicated columns in the
#'  resulting \code{\link{data.frame}}.
#'@param n_iter Number of iterations (default: 5000).
#'@param seed Number for \code{\link{set.seed}}; \code{8} by default. Set to
#'  \code{NULL} for random seed.
#'
#'@details
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
#'the varying values. When so specified, a separate simulation procedure will be
#'run for each combination of the given factors (or, if only one factor is
#'given, for each element of that factor). The
#'\code{\link[POSSA:pow]{POSSA::pow}} function will be able to automatically
#'detect (by default) the factors generated this way in the present
#'\code{\link[POSSA:sim]{POSSA::sim}} function, in order to calculate power
#'separately for each factor combination.
#'
#'@return Returns a \code{\link{data.frame}} that includes the columns
#'  \code{iter} (the iterations of the simulation procedure numbered from
#'  \code{1} to \code{n_iter}), \code{look} (the interim "looks" numbered from
#'  \code{1} to the maximum number of looks, including the final one), and the
#'  information returned by the \code{fun_test} function for H0 and H1 outcomes
#'  (mainly p values; but also other, optional information, if any) and the
#'  corresponding observation numbers.
#'
#'@note
#'
#'For the replicability (despite the randomization), \code{\link{set.seed}} is
#'executed in the beginning of this function, each time it is called; see the
#'\code{seed} parameter.
#'
#' @seealso \code{\link{pow}}
#' @examples
#'
#' # user-defined function to specify sample(s)
#' custom_sample_v1 = function(v1, v2_h) {
#'     list(v1 = rnorm(v1, mean = 0, sd = 10),
#'          v2_h0 = rnorm(v2_h, mean = 0, sd = 10),
#'          v2_h1 = rnorm(v2_h, mean = 5, sd = 10))
#' }
#'
#'
#' # user-defined function to specify significance test(s)
#' custom_test_v1 = function(sampl) {
#'     c(
#'         p_h0 = t.test(sampl$v2_h0, sampl$v1, var.equal = TRUE)$p.value,
#'         p_h1 = t.test(sampl$v2_h1, sampl$v1, var.equal = TRUE)$p.value
#'     )
#' }
#'
#' df_ps_v1 = sim(
#'     fun_obs = custom_sample_v1,
#'     n_obs = c(33, 65, 98),
#'     fun_test = custom_test_v1
#' )
#'
#' @export
sim = function(fun_obs,
               n_obs,
               fun_test,
               n_iter = 5000,
               seed = 8) {
    validate_args(match.call(),
                  list(
                      val_arg(fun_obs, c('function', 'list')),
                      val_arg(fun_test, c('function')),
                      val_arg(n_iter, c('num'), 1)
                  ))

    set.seed(seed)
    # sanity checks for given values
    if (!is.function(fun_obs)) {
        f_s_args = fun_obs
        if (!is.function(f_s_args[[1]])) {
            stop('When "fun_obs" is given as list,',
                 ' the first argument must be a function.')
        }
        fun_obs = f_s_args[[1]]
        f_s_args[[1]] = NULL
        # get all fun_obs argument combinations (different sample types)
        df_combs = sapply(expand.grid(f_s_args), as.vector)
        f_s_a_list = list()
        for (rownum in 1:nrow(df_combs)) {
            f_s_a_list[[rownum]] = as.list(df_combs[rownum,])
        }
    } else {
        # set to have no combinations; single sample test (hence 1 cycle below)
        f_s_a_list = NA
        f_s_args = c()
    }
    if (is.atomic(n_obs)) {
        # if just a vector given, all samples have this as samples sizes
        n_obs_orig = n_obs
        n_obs = list()
        fpars = methods::formalArgs(fun_obs) # required sample names guessed from fun_obs arguments
        for (n_name in fpars[!fpars %in% names(f_s_args)]) {
            n_obs[[n_name]] = n_obs_orig # assign vector to each sample type
        }
    }
    n_look = length(n_obs[[1]])
    if (!all(sapply(n_obs, length) == n_look)) {
        stop('The lengths of the "n_obs" values are unequal.')
    }
    obs_per_it = do.call(Map, c(f = list, n_obs)) # transpose input for iteration below
    obs_names = names(obs_per_it[[1]]) # the names of all sample columns
    list_vals = list() # all end results will be collected in this list
    # set progress bar
    pb = utils::txtProgressBar(
        min = 0,
        max = n_iter * length(f_s_a_list),
        initial = 0,
        style = 3
    )
    pb_count = 0
    for (f_s_a in f_s_a_list) {
        if (is.na(f_s_a[1])) {
            f_s_a = NULL
        }
        for (i in 1:n_iter) {
            pb_count = pb_count + 1
            utils::setTxtProgressBar(pb, pb_count)
            samples = do.call(fun_obs, c(obs_per_it[[n_look]], f_s_a))
            list_vals[[length(list_vals) + 1]] =
                c(
                    iter = i,
                    look = n_look,
                    unlist(f_s_a),
                    unlist(obs_per_it[[n_look]]),
                    fun_test(samples)
                )
            for (lk in (n_look - 1):1) {
                seed = .Random.seed
                for (samp_n in obs_names) {
                    if (endsWith(samp_n, '_h')) {
                        for (h_num in c('0', '1')) {
                            .Random.seed = seed
                            samples[[paste0(samp_n, h_num)]] = sample(samples[[paste0(samp_n, h_num)]],
                                                                      obs_per_it[[lk]][[samp_n]])
                        }
                    } else {
                        .Random.seed = seed
                        samples[[samp_n]] = sample(samples[[samp_n]],
                                                   obs_per_it[[lk]][[samp_n]])
                    }
                }
                list_vals[[length(list_vals) + 1]] =
                    c(
                        iter = i,
                        look = lk,
                        unlist(f_s_a),
                        unlist(obs_per_it[[lk]]),
                        fun_test(samples)
                    )
            }
        }
    }
    close(pb)
    df_pvals = as.data.frame(do.call(rbind, list_vals))
    df_pvals = df_pvals[order(df_pvals$iter, df_pvals$look),]
    for (c_nam in names(n_obs)) {
        class(df_pvals[[c_nam]]) = c(class(df_pvals[[c_nam]]), "possa_n")
    }
    for (fc_nam in names(f_s_args)) {
        class(df_pvals[[fc_nam]]) = c(class(df_pvals[[fc_nam]]), "possa_fac")
    }
    class(df_pvals) = c(class(df_pvals), "possa_df")
    return(df_pvals)
}