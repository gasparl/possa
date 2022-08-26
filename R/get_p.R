#'@title Extract p-value from tests
#'
#'@description This function attempts to extract p values from certain tests
#'  where it could otherwise be complicated to do so. Please make sure, in case
#'  of each new test, whether the function actually returns the values you want.
#'
#'@param x Test results.
#'
#'@details
#'
#' Supported functions: all tests that return the p value as \code{p.value}
#' (including most R \code{stats} test functions, having \code{htest} class),
#' and the \code{aov()} function (\code{aov} and \code{aovlist} classes).
#'
#'@return Returns either a single p value or, in case of multiple p values, a
#'  list or nested list with each p value.
#'
#' @examples
#'
#' get_p(t.test(extra ~ group, data = sleep))
#' # returns 0.07939414
#' # same as printed via t.test(extra ~ group, data = sleep)
#'
#' get_p(prop.test(c(83, 90, 129, 70), c(86, 93, 136, 82)))
#' # returns 0.005585477,
#' # same as printed prop.test(c(83, 90, 129, 70), c(86, 93, 136, 82))
#'
#' get_p(aov(yield ~ block + N * P * K, npk))
#' # returns list of p values
#' # corresponds to summary(aov(yield ~ block + N * P * K, npk))
#'
#' get_p(aov(yield ~ N * P * K + Error(block), npk))
#' # returns nested list of p values (effects per error term)
#' # again corresponds printed p values via summary()
#'
#' @export
get_p = function(x) {
    UseMethod("get_p")
}

#' @export
get_p.default <- function(x) {
    if (is.numeric(x$p.value)) {
        return(x$p.value)
    } else {
        warning(
            'No "p.value" found by the "get_p" function.',
            ' Please double-check your input',
            ' or extract the p value in another way.'
        )
    }
}

#' @describeIn get_p get_p method for class 'aov'
#' @export
get_p.aov <- function(x) {
    terms <- summary(x)[[1]]
    # Trim spaces from the names of the terms
    p = stats::setNames(as.list(terms$`Pr(>F)`), gsub(" ", "", rownames(terms)))
    return(p[-length(p)])
}

#' @describeIn get_p get_p method for class 'aovlist'
#' @export
get_p.aovlist <- function(x) {
    aov_sum = summary(x)
    # Loop over the error strata
    p = list()
    for (i in 1:length(names(aov_sum))) {
        # Get term statistics of the current error stratum
        terms <- aov_sum[[i]][[1]]
        p[[gsub(" ", "", names(aov_sum)[i])]] = stats::setNames(as.list(terms$`Pr(>F)`), gsub(" ", "", rownames(terms)))[-nrow(terms)]
    }
    return(p)
}