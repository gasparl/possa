### Misc internal functions ###

# rounding numbers, for nicer output
ro = function(num,
              round_to = 2,
              leading_zero = FALSE,
              signi = FALSE) {
    validate_args(match.call(),
                  list(val_arg(leading_zero, c('bool'), 1),
                       val_arg(signi, c('bool'), 1)))
    if (is.numeric(num)) {
        value = num
    } else {
        value = as.numeric(as.character(num))
    }
    if (signi == TRUE) {
        formtd = base::formatC(value, format = "fg", digits = round_to)
    } else {
        formtd = format(round(value, round_to),
                        nsmall = round_to,
                        scientific = FALSE)
    }
    formtd = gsub(" ", "", formtd)
    if (leading_zero == FALSE) {
        formtd = sub("0.", ".", formtd, fixed = TRUE)
        tinies = gsub("0", "", formtd, fixed = TRUE) == '.'
        formtd[tinies] = sub('.$', '1',
                             sub(".", "< .", formtd[tinies], fixed = TRUE))
        formtd[value == 1] = "1"
    }
    formtd[value == 0] = "0"
    return(formtd)
}

# check if column exists in the dataframe (when it should)
checkcol = function(df_names, thecols) {
    cols_notfound = c()
    for (colname in thecols) {
        if (!colname %in% df_names) {
            cols_notfound = c(cols_notfound, colname)
        }
    }
    if (length(cols_notfound) > 0) {
        if (length(cols_notfound) ==  1) {
            stop(
                'The column "',
                cols_notfound,
                '" was not found in the data frame. Perhaps check for spelling mistakes.'
            )
        } else {
            stop(
                'The following columns were not found in the data frame: "',
                paste(cols_notfound,
                      collapse = '", "'),
                '". Perhaps check for spelling mistakes.'
            )
        }
    }
}

# check if column exists in the dataframe (when it shouldn't)
name_taken = function(name, dat) {
    if (name %in%  names(dat)) {
        stop(
            'Sorry, the name "',
            name,
            '" is reserved for this function. Remove or rename that column.'
        )
    }
}

# checking suffix consistency
suffixes = function(thenames,
                    feedf) {
    pairs = c()
    for (a_nam in thenames) {
        if (endsWith(a_nam, '_h0')) {
            n_nam = substr(a_nam, 1, nchar(a_nam) - 1)
            if (paste0(n_nam, '1') %in% thenames) {
                pairs = c(pairs, n_nam)
            } else {
                feedf('Found "h0" suffix without matching "h1": ',
                      a_nam)
            }
        } else if (endsWith(a_nam, '_h1')) {
            n_nam = substr(a_nam, 1, nchar(a_nam) - 1)
            if (!paste0(n_nam, '0') %in% thenames) {
                feedf('Found "h1" suffix without matching "h0": ',
                      a_nam)
            }
        }
    }
    return(pairs)
}

## parameter argument valudations

validate_args = function(func_used, evaled_args) {
    feedback = ''
    for (part_feed in evaled_args) {
        feedback = paste0(feedback, part_feed)
    }
    if (feedback != '') {
        func_used = gsub("\\s+", " ", paste(deparse(func_used), collapse = " "))
        feedback = paste0(
            "Arguments are not correct in the '",
            func_used,
            "' function:",
            feedback,
            '\n... Hint: enter help(',
            strsplit(func_used, "\\(")[[1]][1],
            ') for detailed function info.'
        )
        stop(feedback, call. = FALSE)
    }
}

val_arg = function(arg_val,
                   req_types,
                   req_length = 99,
                   # 0 means multiple, 1 means single, all else passes
                   opts = NULL) {
    failed = FALSE
    arg_name = paste(deparse(substitute(arg_val)), collapse = "")
    if (length(arg_val) > 1) {
        if (req_length == 1 &&
            !is.list(arg_val) && !is.data.frame(arg_val)) {
            failed = TRUE
        }
    } else if (req_length == 0) {
        failed = TRUE
    }
    valid_types = c('char', 'num', 'bool', 'null', 'df', 'list', 'function')
    if (!all(req_types %in% valid_types)) {
        stop(
            'invalid req_types: ',
            paste(req_types, collapse = ', '),
            '\nshould be: ',
            paste(valid_types, collapse = ', ')
        )
    }
    req_types = replace(req_types, req_types == 'char', 'character')
    req_types = replace(req_types, req_types == 'num', 'double')
    req_types = replace(req_types, req_types == 'bool', 'logical')
    req_types = replace(req_types, req_types == 'null', 'NULL')
    req_types = replace(req_types, req_types == 'df', 'data.frame')
    if ((!typeof(arg_val) %in% req_types)
        && (!('data.frame' %in% req_types &&
              is.data.frame(arg_val)))
        && (!('function' %in% req_types &&
              is.function(arg_val)))
        && (!('double' %in% req_types &&
              is.integer(arg_val))) &&
        (!('list' %in% req_types &&
           is.list(arg_val)))) {
        failed = TRUE
    } else if (typeof(arg_val) == 'character' &&
               (!is.null(opts)) && (!arg_val %in% opts)) {
        failed = TRUE
    }
    if (failed == TRUE) {
        req_types = replace(req_types,
                            req_types == 'character',
                            '"character" (string)')
        req_types = replace(req_types, req_types == 'double', '"double" (numeric)')
        req_types = replace(req_types,
                            req_types == 'logical',
                            '"logical" (boolean)')
        req_types = replace(req_types, req_types == 'data.frame', '"data.frame"')
        if (!is.null(opts)) {
            if (suppressWarnings(all(!is.na(as.numeric(opts))))) {
                opts_add = paste0(
                    ' The only acceptable strings or numbers are "',
                    paste(opts, collapse = '", or "'),
                    '".'
                )
            } else {
                opts_add = paste0(
                    ' The only acceptable strings are "',
                    paste(opts, collapse = '", or "'),
                    '".'
                )
            }
        } else {
            opts_add = ''
        }
        if (req_length == 1) {
            len_add = ' must not be a vector, and'
        } else if (req_length == 0) {
            len_add = ' must be a vector, and'
        } else {
            len_add = ''
        }
        arg_msg = paste0(
            '\nThe argument "',
            arg_name,
            '"',
            len_add,
            ' must be ',
            paste(req_types, collapse = ', or '),
            '.',
            opts_add
        )
        return(arg_msg)
    } else {
        return('')
    }
}
