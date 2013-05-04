# Copyright: Marc Wildi
# 15.01.2012
# http://blog.zhaw.ch/idp/sefblog

#_____________________________________________________________


#gdp_q<-gdp
#x<-xh


# New 2012-code: computes spectral estimates based on DFT
spec_comp<-function(insamp,x,d)
{
#insamp<-260-anf+1
#insamp<-len
#  K<-(insamp-1)/2
  weight_func<-NULL#matrix(nrow=K+1,ncol=dim(x)[2])

  if (d==1)
  {
    K<-length(periodogram_bp(diff(x[1:insamp,1]), 1, insamp-1)$fourtrans)-1
    weight_func<-1:(K+1)
#    if (!GDP_T)
#    {
      weight_func<-periodogram_bp(diff(x[1:insamp,1]), 1, insamp-1)$fourtrans
#    } else
#    {
#      weight_func<-periodogram_bp(diff(gdp_q[1:round(insamp/3)]), 1, round(insamp/3)-1)$fourtrans
#      weight_func<-c(weight_func,rep(0,K+1-length(weight_func)))
#      weight_func<-weight_func*exp(-1.i*(0:K)*pi*publication_lag/K)
#    }
    # explaining variables                            ts.plot(x[1:insamp,1])
    if (length(weight_func)>1)
    {
      for (j in 2:length(x[1,]))  #j<-2
      {
    # Since the data is integrated one uses the pseudo-periodogram: diff(data) and d=1
        weight_func<-cbind(weight_func,periodogram_bp(diff(x[1:insamp,j]), 1, insamp-1)$fourtrans)
      }
    }
  } else
  {
    weight_func<-periodogram_bp((x[1:insamp,1]), 0, insamp)$fourtrans
    K<-length(periodogram_bp(x[1:insamp,1], 0, insamp)$fourtrans)-1
#    if (!GDP_T)
#    {
      weight_func<-periodogram_bp(x[1:insamp,1], 0, insamp)$fourtrans
#    } else
#    {
#      weight_func<-periodogram_bp(gdp_q[1:round(insamp/3)], 0, round(insamp/3))$fourtrans
#      weight_func<-c(weight_func,rep(0,K+1-length(weight_func)))
#      weight_func<-weight_func*exp(-1.i*(0:K)*pi*publication_lag/K)
#    }


    # explaining variables                            ts.plot(x[1:insamp,1])
    if (length(weight_func)>1)
    {
      for (j in 2:length(x[1,]))  #j<-2
      {
    # Since the data is integrated one uses the pseudo-periodogram: diff(data) and d=1
        weight_func<-cbind(weight_func,periodogram_bp((x[1:insamp,j]), 0, insamp)$fourtrans)
      }
    }

  }
  dimnames(weight_func)[[2]]<-dimnames(x)[[2]]
  #weight_func[,1]<-periodogram_bp(diff(gdp[1:insamp]), 1, insamp-1)$fourtrans
  # if i1<-T then weight_constraint imposes corresponding values of amplitude functions in frequency zero

#  ts.plot(abs(weight_func)[,1])
  return(list(weight_func=weight_func))
}




# DFT (old code but still in use for new 2012-version...)
periodogram_bp <- function(x, dd, n.pg)
{
    ## Preparations
    n.fit  <- length(x)
    xx     <- x[((n.fit-n.pg+1):n.fit)]
    npg2   <- (n.pg/2)
    perall <- 0*0:npg2
    fourtrans<-perall
    
    ## Case without a seasonal component
    if (dd < 3)
    {
        for (j in (1:npg2)) #j<-1
        {
            fourtrans[j+1] <- xx%*%exp((1:(2*npg2))*1.i*j*pi/npg2)
            term2 <- (1-exp(j*1.i*pi/npg2))^(dd)
            fourtrans[j+1] <- fourtrans[j+1]/(1-min(dd,1)+min(1,dd)*term2)
            
            perall[j+1] <- abs(fourtrans[j+1])^2
        }
    }
    
    ## Case with a seasonal component, special treatment for Pi/6
    if (dd >= 3)
    {
        for (j in (1:npg2)[(-npg2/6)*(1:6)])
        {
            fourtrans[j+1] <- xx%*%exp((1:(2*npg2))*1.i*j*pi/npg2)
            term2 <- abs(1-exp(j*1.i*pi/npg2))^2
            term3 <- abs(1-exp(12*j*1.i*pi/npg2))^2
            perall[j+1] <- abs(fourtrans[j+1])/(term2*term3)
        }
        perall[(npg2/6)*(1:6)+1] <- max(perall)*100000
    }
    
    ## Output
    return(list(perall=perall,fourtrans=fourtrans))
}




## DFT (old code but still in use for new 2012-version...)
#periodogram_bp <- function(x, dd, n.pg) {
#    n.fit  <- length(x)
#    xx     <- tail(x, n.pg)
#    npg2   <- n.pg / 2
#    
#    ## Case with dd == 0
#    if (dd == 0) {
#        M <- exp(outer(1:n.pg, 1:npg2) * 1.i * pi / npg2)
#        fourtrans <- c(0, xx %*% M)
#        perall <- abs(fourtrans)^2
#    }
#
#    ## Case without a seasonal component
#    if (dd > 0 & dd < 3) {
#        M <- exp(outer(1:n.pg, 1:npg2) * 1.i * pi / npg2)
#        fourtrans <- c(0, xx %*% M)
#        for (j in 1:npg2) {
#            fourtrans[j+1] <- xx %*% exp((1:n.pg)*1.i*j*pi/npg2)
#            term2 <- (1-exp(j*1.i*pi/npg2))^(dd)
#            fourtrans[j+1] <- fourtrans[j+1]/(1-min(dd,1)+min(1,dd)*term2)
#            perall[j+1] <- abs(fourtrans[j+1])^2
#        }
#    }
#
#    ## Case with a seasonal component, special treatment for Pi/6
#    if (dd >= 3) {
#        perall <- rep(0, npg2 + 1)
#        fourtrans<-perall
#        for (j in (1:npg2)[(-npg2/6)*(1:6)])
#          {
#            fourtrans[j+1] <- xx%*%exp((1:(2*npg2))*1.i*j*pi/npg2)
#            term2 <- abs(1-exp(j*1.i*pi/npg2))^2
#            term3 <- abs(1-exp(12*j*1.i*pi/npg2))^2
#            perall[j+1] <- abs(fourtrans[j+1])/(term2*term3)
#          }
#        perall[(npg2/6)*(1:6)+1] <- max(perall)*100000
#      }
#
#    ## Output
#    return(list(perall = perall, fourtrans = fourtrans))
#  }
