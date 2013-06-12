#' Create a matrix appropriate for passing to mdfa().
#' 
#' Creates a matrix suitable for use in an imdfa model fit.
#' The formula can be structured like this "y ~ y + x" or in the more standard
#' "y ~ x" format, as desired.
#' 
#' @param formula a formula describing the desired model
#' @param a matrix, dataframe, list or xts object
#' @export
mdfa_model_matrix <- function(formula, data) {
    tm <- terms(formula, data=data)
    vnames <- dimnames(attr(tm, "factors"))
    fm <- reformulate(c(vnames[[1]][1],vnames[[2]]), response = NULL, intercept = FALSE)
    mat <- model.matrix(fm, data=as.data.frame(data, stringsAsFactors = FALSE))
    mat <- cbind(mat[,which(vnames[[1]][1] == vnames[[2]])], mat)
    colnames(mat) <- c(vnames[[1]][1],vnames[[2]])
    return(mat)
}

#' Fit an i-mdfa model.
#'
#' Main function to calculate filter using I-MDFA code. If a formula is not used, the first column of
#' the data is the dependent, and the remaining columns are the predictors.
#'
#' @param data a data.frame, list, matrix, vector or xts data series
#' @param formula a symbolic description of the model as a formula object. If a formula is used,
#' the data object must be a data frome or list.
#' @param ... further arguments to be passed to mdfa_core
#' @keywords dfa mdfa imdfa
#' @export
#' @examples
#' m <- mdfa(y ~ a + b, df)
#' summary(m)
#' plot(m)
#' coef(m)
#' mp <- predict(m)
#' plot(mp)
mdfa.default <- function(data, formula = NULL, keep_data = TRUE, spectral_estimate = NULL, d = 0,
                         cutoff = pi/12, Gamma = NULL, L = 24, ...) {
    args_list <- as.list(environment())
    cl <- match.call()
    tm <- NULL
	cls_data = 't'
    if (is.xts(data)) {
        ix_data <- index(data)
        data <- coredata(data)
        cls_data <- 'xts'
    }
    if (!is.null(formula)) {
        data <- mdfa_model_matrix(formula, data)
    }
	if (!is.matrix(data)){
        if (is.list(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
        data <- as.matrix(data)
    }
    if (NCOL(data) == 1) data <- cbind(data,data) 
    if (is.null(spectral_estimate)) {
        spectral_estimate <- calc_dfts(data, NROW(data), d)
        args_list[['spectral_estimate']] <- spectral_estimate
    }
    mdfa.orig <- do.call(mdfa_core, args_list)
    mdfa.orig$call <- cl
    mdfa.orig$formula <- formula
    if (cls_data == 'xts') data <- xts(data, order.by = ix_data)
    if (keep_data) mdfa.orig$data <- data
    class(mdfa.orig) <- c('mdfa', class(mdfa.orig))
    if (!is.null(formula)){
        colnames(mdfa.orig$b) <- attr(terms(formula),'term.labels')
    } else {
        colnames(mdfa.orig$b) <- colnames(data)[-1]
    }
    class(mdfa.orig$b) <- c('mdfa_coef', class(mdfa.orig$b))
    return(mdfa.orig)
}

#' Print methods for class "mdfa" and "mdfa_coef".
#'
#' Prints out the str() representation of the mdfa object. For the `mdfa_coef`
#' object, prints the matrix of coefficients.
#'
#' @param object fitted mdfa object
#' @S3method print mdfa
print.mdfa <- function(x, ...) {
	invisible(str(x))
}
print.mdfa_coef <- function(x, ...) {
	NextMethod(x)
}

#' Core calculation for mdfa
#' 
#' Wrapper around Marc's code that sets defaults and is available for dispatch 
#' from within the mdfa family of functions.
#' 
mdfa_core <- function(L, Gamma, cutoff, spectral_estimate, K=NROW(spectral_estimate)-1, lambda=0, Lag=0, expweight=0, 
                      i1=FALSE, i2=FALSE, weight_constraint=1, lambda_cross=0, lambda_decay=c(0, 0), lambda_smooth=0.1,
                      lin_expweight=FALSE, shift_constraint=0, grand_mean=TRUE, ...) {
    args_list <- as.list(environment())
    args_list[['data']] <- NULL
    args_list[['weight_func']] <- spectral_estimate
    args_list[['spectral_estimate']] <- NULL
    if (is.null(Gamma))
        args_list[['Gamma']] <- lowpass_filter_spec(K, L, cutoff)
    return(do.call(mdfa_analytic_new, args_list))
}


mdfa <- function(x, ...) {
    UseMethod("mdfa", x)
}

mdfa.formula <- function(fmla, data, keep_data = TRUE, d = 0, ...) {
    args_list <- as.list(environment())
    cl <- match.call()
    tm <- NULL
    cls_data = 't'
    if (is.xts(data)) {
        ix_data <- index(data)
        data <- coredata(data)
        cls_data <- 'xts'
    }
    data <- mdfa_model_matrix(formula, data)
    if (!is.matrix(data)){
        if (is.list(data)) data <- as.data.frame(data, stringsAsFactors = FALSE)
        data <- as.matrix(data)
    }
    if (NCOL(data) == 1) data <- cbind(data,data) 
    spectral_estimate <- calc_dfts(data, NROW(data), d)
    args_list[['spectral_estimate']] <- spectral_estimate
    mdfa.orig <- do.call(mdfa_core, args_list)
    mdfa.orig$call <- cl
    mdfa.orig$formula <- formula
    if (cls_data == 'xts') data <- xts(data, order.by = ix_data)
    if (keep_data) mdfa.orig$data <- data
    class(mdfa.orig) <- c('mdfa', class(mdfa.orig))
    colnames(mdfa.orig$b) <- attr(terms(formula),'term.labels')
    class(mdfa.orig$b) <- c('mdfa_coef', class(mdfa.orig$b))
    return(mdfa.orig)
}

mdfa.data.frame <- function(data, keep_data = TRUE, d = 0, ...) {
    args_list <- as.list(environment())
    cl <- match.call()
    data <- as.matrix(data)
    if (NCOL(data) == 1) data <- cbind(data,data)
    spectral_estimate <- calc_dfts(data, NROW(data), d)
    args_list[['spectral_estimate']] <- spectral_estimate
    mdfa.orig <- do.call(mdfa_core, args_list)
    mdfa.orig$call <- cl
    if (keep_data) mdfa.orig$data <- data
    class(mdfa.orig) <- c('mdfa', class(mdfa.orig))
    class(mdfa.orig$b) <- c('mdfa_coef', class(mdfa.orig$b))
    return(mdfa.orig)
}

mdfa.xts <- function(data, keep_data = TRUE, d = 0, ...) {
    warning("mdfa.xts not yet implemented")
}

mdfa.matrix <- function(data, keep_data = TRUE, d = 0, ...) {
    args_list <- as.list(environment())
    cl <- match.call()
    if (NCOL(data) == 1) data <- cbind(data,data)
    spectral_estimate <- calc_dfts(data, NROW(data), d)
    args_list[['spectral_estimate']] <- spectral_estimate
    mdfa.orig <- do.call(mdfa_core, args_list)
    mdfa.orig$call <- cl
    if (keep_data) mdfa.orig$data <- data
    class(mdfa.orig) <- c('mdfa', class(mdfa.orig))
    class(mdfa.orig$b) <- c('mdfa_coef', class(mdfa.orig$b))
    return(mdfa.orig)
}

mdfa.spectral_estimate <- function(spectral_estimate, L, Gamma, cutoff, ...) {
    args_list <- as.list(environment())
    cl <- match.call()
    mdfa.orig <- do.call(mdfa_core, args_list)
    mdfa.orig$call <- cl
    class(mdfa.orig) <- c('mdfa', class(mdfa.orig))
    class(mdfa.orig$b) <- c('mdfa_coef', class(mdfa.orig$b))
    return(mdfa.orig)
}
