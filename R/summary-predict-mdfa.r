#' Summary method for class "predict_mdfa".
#'
#' Long description
#'
#' @param object mdfa prediction object
#' @S3method summary predict_mdfa
summary.predict_mdfa <- function(object) {
    warning("Summary not yet implemented for predict_mdfa.")
}

#' Print method for class "predict_mdfa".
#'
#' Long description
#'
#' @param object mdfa prediction object
#' @S3method print predict_mdfa
print.predict_mdfa <- function(object) {
	NextMethod(print, object)
}