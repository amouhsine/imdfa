#' Diagnostic plots for an mdfa model.
#'
#' Plot amplitude and phase response on the transfer function for a fitted MDFA model
#'
#' @param object fitted mdfa object
#' @method plot mdfa
#' @S3method plot mdfa
#' @export
plot.mdfa <- function(object) {
    trffkt <- object$trffkt
    nc <- NCOL(trffkt)
    K <- object$options$K
    
    nms <- dimnames(object$data)[[2]]
    
    ampl <- NULL
    phas <- NULL
    for (i in 1:nc) {
        ampl <- cbind(ampl, abs(trffkt[,i])^2)
        phas <- cbind(phas, -Arg(trffkt)[2:(K+1),i]/((1:(K))*pi/K))
    }

    # label charts from 0:pi on 1/6
    #mw uses: axis(1,at=1+0:6*N/12,labels=c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi"))
    labels_at <- 1 + 0:6*K/6
    lbl_text <- c("0","pi/6","2pi/6","3pi/6","4pi/6","5pi/6","pi")
    
    par(mfrow=c(2,1))
    ##amplitude function
    plot.zoo(ampl, plot.type='single', col=1:nc, xlab='Frequency', ylab='Amplitude', xaxt='n')
    axis(side=1, at=labels_at, labels=lbl_text)
    
    ##TODO: need to pull names from the formula to get the legend to work properly
    #legend('topright', legend=nms, text.col=1:nc, ncol=2)
    
    ##phase function
    ##mw uses: plot(-Arg(i_mdfa_opt$trffkt)[2:(K+1),1]/((1:(K))*pi/K),type="l",axes=F,col="blue")
    plot.zoo(phas, plot.type='single', col=1:nc, xlab='Frequency', ylab='Delay', xaxt='n')
    axis(side=1, at=labels_at, labels=lbl_text)
    
}

#' Plot mdfa model coefficients
#' 
#' Plot the coefficients directly. One line for each dependent variable.
#' 
#' @S3method plot mdfa_coef
#' 
plot.mdfa_coef <- function(object) {
    dat <- melt(object, varnames = c("Time", "Variable"), value.name = "Coefficient")
    ggplot(dat, aes(x = Time, y = Coefficient)) +
            geom_line(aes(color = Variable))
}