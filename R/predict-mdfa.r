#' Make predictions based on a fitted mdfa model.
#'
#' Long description
#'
#' @param object mdfa object
#' @param data new data for use in prediction
#' @export
predict.mdfa <- function(object, data=NULL) {
    if (is.null(data)) data <- object$data[,-1]
    cls_data <- class(data)
    if (is.xts(data)) {
        ix_data <- index(data)
        data <- coredata(data)
        cls_data <- 'xts'
    }
    if (!is.null(object$formula) & !is.null(data)) {
        data <- mdfa_model_matrix(object$formula, data)[,-1]
    }
    nr <- NROW(data)
    nc <- NCOL(data)
    b <- object$b
    if (nc != NCOL(b)) stop("Number of columns in data does not match number of columns in the model.")
    xff <- matrix(NA, ncol=nc, nro=nr)
    for (i in 1:nc){
        xff[,i] <- filter(data[,i], b[,i], method="convolution", sides=1)
    }
    pr <- rowSums(xff)
    if (cls_data == 'xts') {
        pr <- xts(pr, order.by = ix_data)
    }
	class(pr) <- c(class(pr), "predict_mdfa")
	return(pr)
}

#' Print method for class "predict_mdfa".
#'
#' Long description
#'
#' @param object mdfa prediction object
#' @export
print.predict_mdfa <- function(x, ...) {
    invisible(print(as.vector(x)))
}