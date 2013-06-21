## simple wrapper around marc's standard periodogram function
## data_series must be a matrix or matrix-alike (not a vector)
calc_dfts <- function(input_srs, iis_rows=NULL, seasonality=0) {
    nr <- iis_rows
    if (is.null(iis_rows)) nr <- NROW(input_srs)

    weight_func <- NULL
    for (j in 1:NCOL(input_srs)) {
        wc <- periodogram_bp(input_srs[,j], seasonality, nr)$fourtrans
        weight_func <- cbind(weight_func, wc)
    }
    weight_func
}

## straight copy/paste of marc's standard lowpass filter spec
lowpass_filter_spec <- function(grid_size, filter_length, noise_cutoff) {
    ord <- filter_length
    cutoff <- noise_cutoff
    
    b<-c(cutoff/pi,(1/pi)*sin(cutoff*1:ord)/(1:ord))
    b<-b/(sum(b)+sum(b[2:(ord+1)]))
    
    trffkt<-0:grid_size
    for (k in 0:grid_size) {
        omegak<-k*pi/grid_size
        trffkt[k+1]<-b%*%(cos(omegak*0:ord)*c(1,rep(2,ord)))
    }
    
    Gamma <- abs(trffkt)
    return(Gamma)
}

## code formatting for Marc's code.
## Here so that there is a record of how it was performed.
#library(formatR)
#tidy.source(source = 'imdfa/inst/I-MDFA.r',
#        keep.comment = FALSE, keep.blank.line = FALSE,
#        reindent.spaces = 4,
#        file = 'imdfa/R/I-MDFA_tidy.r')


