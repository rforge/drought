#' Multivariate frequency analysis with copula (return period with exceedance probability)
#' Useing Gumbel copula as an example
#' @param X are the drought properties or indices
#' @param Y are the drought properties or indices
#' @param EL is the average reocurrence time
#' @return The multivariate return period
#' @export
#' @examples
#' X=runif(120, min = 0, max = 100)+100
#' Y=runif(120, min = 0, max = 100)+100
#'CopFreq(X,Y,EL=1)
 
CopFreq<-function (X,Y,EL) 
{
  
  
  # Assume the mean interval EL=1 (or with annual maxima)

  # Compute the univariate return period
  
  #F1=stats::ecdf(X)
  #EP1=F1(X);
  #plot(X,EP1,type='s',xlab="X", ylab="Pro.")
  #F2=stats::ecdf(Y)
  #EP2=F2(Y);
  #plot(Y,EP2,type='s',xlab="Y", ylab="Pro.")
  
  pa<-MASS::fitdistr(X,"gamma")
  X0=seq(min(X)*0.8,max(X)*1.2,(max(X)-min(X))/100)
  P1=stats::pgamma(X0,pa$estimate[1],pa$estimate[2])
  T1=1/(1-P1); # Exceedance return period
  
  
  pa<-MASS::fitdistr(Y,"gamma")
  Y0=seq(min(Y)*0.8,max(Y)*1.2,(max(Y)-min(Y))/100)
  P2=stats::pgamma(Y0,pa$estimate[1],pa$estimate[2])
  T2=1/(1-P2);# Exceedance return period
  
  
  # Compute the joint return period
  
  gumbel.cop <- copula::gumbelCopula(4, dim=2)
  
  #u <- copula::rCopula(100, gumbel.cop)
  
  u=copula::pobs(cbind(X,Y))
  
  theta <- copula::fitCopula(gumbel.cop, u, method="itau")
  
  #cop <- copula::onacopula("Gumbel", list(theta,1:2))

  
  uu<-matrix(seq(0.01,0.99,0.01), nrow=99,ncol=99)
  MTA<-matrix(seq(0.01,0.99,0.01), nrow=99,ncol=99)
  MTO<-matrix(seq(0.01,0.99,0.01), nrow=99,ncol=99)
  
  for (i in 1:99)
  {
    for (j in 1:99)
    {
    p=copula::pCopula(uu[i,j],gumbel.cop) 
      
    MTA[i,j]= 1/(1-uu[i,j]-uu[i,j]+p)*EL  
    MTO[i,j]= 1/(1-p)*EL  
    
    }
  }
  
  
  print("The Univariate and multivariate return period:  AND and OR case")
  
  
  
  result<-list(UT1=T1,UT2=T2,MTA=MTA,MTO=MTO)
  
  return(result)
  
}
