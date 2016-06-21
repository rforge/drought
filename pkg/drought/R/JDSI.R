#' Compute the multivariate drought index with joint distribution 
#' 
#' @param X  is The vector of a monthly hydro-climatic variable of n years. 
#' @param Y  is The vector of a monthly hydro-climatic variable of n years.
#' @param ts is the accumulated time scale. 
#' @return The multivariate drought index of different time scales from the marginal probability (or percentile) 
#' @export
#' @examples
#' X=runif(120, min = 0, max = 100)
#' Y=runif(120, min = 0, max = 100)
#' fit<-JDSI(X,Y,ts=6)  
#' fit$JDSI_A 
JDSI<-function (X,Y,ts=6) 
{
  
  X=matrix(X,ncol=1)
  Y=matrix(Y,ncol=1)
  
  # Obtain the accumulation of the variable for a spec`ific time scale (ts) 
  
  AX=ACCU(X,ts)
  AY=ACCU(Y,ts)
  
  #  define the SI=Y0 for each month
  #  define the empirical distribution
  

  
  JDSI_A=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  JDSI_O=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  JDSI_K=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  
  for  (k in 1:12)
  {
    
    mdx=AX[,k]
    mdy=AY[,k]
    
    # if the first number is NA, exclude the first number and then compute SPI
    
    if (k<ts) 
    {
      xd=mdx[2:length(mdx)]
      yd=mdy[2:length(mdy)]
      
      
      p=BiEmp(xd,yd);
      u=copula::pobs(cbind(xd,yd))
      
      
      p1=p
      p2=u[,1]+u[,2]-p
      p3=copula::Kn(p,u)*length(xd)/(length(xd)+1)
      
      
      JDSI_A[2:length(mdx),k]=stats::qnorm(p1)
      JDSI_O[2:length(mdx),k]=stats::qnorm(p2)
      JDSI_K[2:length(mdx),k]=stats::qnorm(p3)
      
    }
    
    else
      # if the first number is not NA, take the whole month to compute SPI
    {  
      
      xd=mdx
      yd=mdy    
      
      p=BiEmp(xd,yd);
      u=copula::pobs(cbind(xd,yd))
      
      p1=p
      p2=u[,1]+u[,2]-p
      p3=copula::Kn(p,u)*length(xd)/(length(xd)+1)

      JDSI_A[,k]=stats::qnorm(p1)
      JDSI_O[,k]=stats::qnorm(p2)
      JDSI_K[,k]=stats::qnorm(p3)
    }
    
  }
  
  
  JDSI1=as.vector(c(t(JDSI_A)))
  JDSI2=as.vector(c(t(JDSI_O)))
  JDSI3=as.vector(c(t(JDSI_K)))

  
  result<-list(JDSI_A=JDSI1,JDSI_O=JDSI2,JDSI_K=JDSI3)
  
  return(result)
}
