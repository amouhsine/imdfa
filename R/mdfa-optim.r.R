#' mdfa set up for optimization
#' 
#' WARNING: This function is very alpha, parameters may change.
#' A wrapper for mdfa that takes a vector of parameters
#' so that it can be used in optim() calls. The order of parameters in x:
#' 1. cutoff
#' 2. lambda
#' 3. expweight
#' 4. lambda_cross
#' 5. lambda_decay[[1]]
#' 6. lambda_decay[[2]]
#' 7. lambda_smooth
#' 
#' If the criteron is a function, it should take a mdfa object
#' and return a numeric vector of length one.
#' 
#' @param x a vector of parameters to optimize on
#' @param spectral_estimate a spectral estimate object
#' @param Gamma
#' @param criterion a string selecting one of the fit statistics or a function returning one number
#' @param crit_args a list of further named arguments to pass to criteron
#' @param ... further arguments to pass to mdfa_core
#' @export
mdfa_optim <- function(x, L, spectral_estimate, Gamma, criterion, crit_args=NULL, ...) {
    args_list <- c(as.list(environment()), list(...))
    args_list$x <- NULL
    args_list$cutoff <- x[1]
    args_list$lambda <- x[2]
    args_list$expweight <- x[3]
    args_list$lambda_cross <- x[4]
    args_list$lambda_decay <- x[5:6]
    args_list$lambda_smooth <- x[7]
    m <- do.call(mdfa_core, args_list)
    if(is.function(criterion)) {
        do.call(criterion, c(m, crit_args))
    } else {
        mdfa_orig[[criterion]]
    }
}
