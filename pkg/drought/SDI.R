#' Compute the standardized drought index (will be replaced as UMDI)
#' 
#' Based on the vector of a monthly hydro-climatic variable, the standardized drought 
#' index is computed.Note here the standardized precipitation index (SPI) is used as 
#' the example of the the drought index in the univaraite case. It also represents other drought indices compuated in the similar way as SPI.
#' 
#' Apart from the standardized drought index, the percentile (probability) is also provided,
#' considering the percentile, such as soil moisture percentile, is also used for drought characterization. 
#' 
#' @param X The vector of a monthly hydro-climatic variable of n years. ts is the accumulated time scale.
#' @return The univariate and multivariate drought index of different time scale from both the empirical and gamma distribution
#' @export
#' @examples
#' X=runif(120, min = 0, max = 100)
#' fit<-SDI(X,ts=3) # Compute the 3 month drought index
#' fit$Probemp # Get the empirical probability
#' fit$SDIemp # Get the empirical drought index 

SDI<-function (X,ts=6,dist="emp") 
{
  
  #X=runif(120, min = 0, max = 100)
  
  #X=c(74.1 ,135.3 , 69.3  ,21.1  ,12.4 , 69.8 , 58.9  ,58.4  ,50.3  ,23.4 , 49.6 ,163.3  ,31.7 , 16.8, 108.7 , 30.8  ,23.7 , 31.1  ,34.2  ,19.0 , 11.9 , 26.5 ,128.8,  79.6  ,30.4 , 35.9 , 45.0  , 7.4  ,23.5 ,109.0 , 23.8  ,87.9 , 22.8  ,54.7  ,39.1 , 80.7  ,36.2 , 92.6,  22.9 ,170.7 , 91.8 ,396.1 ,28.1 ,5.8 ,22.9 ,64.2,70.1 ,85.9)
  
  ts=6
  
  X=matrix(X,ncol=1)
  
  if ((length(X)/12 != round(length(X)/12)))
  {
    return("ERROR: The input data should be the vector of the monthly data for a few years")
    
  }
  
  
  if (dist!= "emp" & dist != "gam") 
  {
    stop("Only the empirical and gamma distributions are included for the current version! 
         Other distribution form  will be updated soon.
         Please select either emp or gam")
  }
  
  # Obtain the accumulation of the variable for a specific time scale (ts) 
  
  XA=ACCU(X,ts)
  
  #  define the SI=Y0 for each month
  #  define the empirical distribution
  
  P_emp=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  P_gam=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  SI_emp=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  SI_gam=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  
  for  (k in 1:12)
  {
    md=XA[,k]
    
    # if the first number is NA, exclude the first number and then compute SPI
    
    if (is.na(md[1])==TRUE) 
    {
      xd=md[2:length(md)]
      
      # Get the empirical probability (Weibull)
      empcdf <- ecdf(xd) 
      cdf0 <- sapply(xd, empcdf)
      cdf0<-cdf0*length(xd)/(length(xd)+1)
      # Get the gamma probability
       par<-MASS::fitdistr(xd,"gamma")
      #par<-fitdistr(xd,"gamma")
      gam_cdf=pgamma(xd,par$estimate[1],par$estimate[2])
      
      #  Save the probability and obtain the index
      
      P_emp[2:length(md),k]=matrix(cdf0,ncol=1)
      P_gam[2:length(md),k]=matrix(gam_cdf,ncol=1)
      
      SI0 <- qnorm(cdf0)
      SI1 <- qnorm(gam_cdf)
      
      SI_emp[2:length(md),k]=matrix( SI0,ncol=1)
      SI_gam[2:length(md),k]=matrix(SI1,ncol=1)
      
    }
    
    else
      # if the first number is not NA, take the whole month to compute SPI
    {  xd=md
       
       # Get the empirical probability (Weibull)
       empcdf <- ecdf(xd) 
       cdf0 <- sapply(xd, empcdf)
       cdf0<-cdf0*length(xd)/(length(xd)+1)
       
       # Get the gamma probability
       par<-MASS::fitdistr(xd,"gamma")
       
       gam_cdf=pgamma(xd,par$estimate[1],par$estimate[2])
       
       #  Save the probability and obtain the index
       
       P_emp[,k]=matrix(cdf0,ncol=1)
       P_gam[,k]=matrix(gam_cdf,ncol=1)
       
       SI0 <- qnorm(cdf0)
       SI1 <- qnorm(gam_cdf)
       
       SI_emp[,k]=matrix( SI0,ncol=1)
       SI_gam[,k]=matrix(SI1,ncol=1)
    }
    
  }
  
  # result<-list(Probemp=P_emp,SDIemp=SI_emp,Probgam=P_gam,SDIgam=SI_gam)
  if (dist=="emp") 
    
  {result<-list(Prob=P_emp,SDI=SI_emp) }
  
  else if (dist =="gam")
    
  {result<-list(Prob=P_gam,SDI=SI_gam) }
  
  else
  {
    stop("Only the empirical and gamma distributions are included! Please select either emp or gam ")
  }
  
  return(result)
}
