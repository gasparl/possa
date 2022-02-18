#'@title Simulation procedure
#'
#'@description This function performs the simulation procedure in order to get p
#'  values. The sample is simulated via the given \code{f_sample} function, and
#'  the significance test is performed via the given \code{f_test} function. The
#'  numbers of observations per look (for a sequential design) are specified in
#'  \code{n_obs}.
#'@param f_sample Function. create any sample variables; use list for varying factors and indicate different parameters in the list elements (as number vectors)
#'@param n_obs A number, a numeric vector, or a list of numeric vectors. (Any of the numbers can always be \code{NA} values as well.)
#' Specifies samples sizes for any sample variable from f_sample that is to be used in f_test. H1 and H2 has to be specified with "_h" ending, then it's converted to "_h0" and "_h1".
#'@param f_test Function. Uses as input the sample variables created in f_sample AND specified in n_obs. The obtained p values must be specified with "p_" start and "_h0" and "_h1" endings (e.g., p_h0, p_h1; p_ttest_h0, p_ttest_h1).
#'@param n_iter Number of iterations.
#'
#'@details
#'
# For details about blah.
#'
#'@return Returns \code{\link{data.frame}}, etc.
#'
#'@note
#'
#' Some notes.
#'
#'@references
#'
#'Lakens, D. (2014). Performing high-powered studies efficiently with sequential
#'analyses: Sequential analyses. European Journal of Social Psychology, 44(7),
#'701–710. \doi{https://doi.org/10.1002/ejsp.2023}
#'
#' @seealso \code{\link{pow}}
#' @examples
#' # some sim
#'
#' @export
sim = function(f_sample,
               n_obs,
               f_test,
               n_iter = 1000,
               seed = 8) {
    validate_args(match.call(),
                  list(
                      val_arg(f_sample, c('function', 'list')),
                      val_arg(f_test, c('function')),
                      val_arg(n_iter, c('num'), 1)
                  ))

    set.seed(seed)
    # sanity checks for given values
    if (!is.function(f_sample)) {
        f_s_args = f_sample
        if (!is.function(f_s_args[[1]])) {
            stop('When "f_sample" is given as list,',
                 ' the first argument must be a function.')
        }
        f_sample = f_s_args[[1]]
        f_s_args[[1]] = NULL
        # get all f_sample argument combinations (different sample types)
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
        fpars = formalArgs(f_sample) # required sample names guessed from f_sample arguments
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
            setTxtProgressBar(pb, pb_count)
            samples = do.call(f_sample, c(obs_per_it[[n_look]], f_s_a))
            list_vals[[length(list_vals) + 1]] =
                c(
                    iter = i,
                    look = n_look,
                    unlist(f_s_a),
                    unlist(obs_per_it[[n_look]]),
                    f_test(samples)
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
                        f_test(samples)
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