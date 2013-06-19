#' Plot mdfa predictions.
#'
#' Long description
#'
#' @param object an mdfa prediction object
#' @param window time period of plot
#' @S3method plot predict_mdfa
plot.predict_mdfa <- function(object, orig_data=NULL, window="all") {
	if (is.null(orig_data)) {
        
    } else if (is.xts(orig_data)) {
        
    } else {
        
    }
}