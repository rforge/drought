#' Compute the multivariate drought index
#' 
#' Based on the vector of a monthly hydro-climatic variable, 
#' the multivariate drought index is compuated based on the margional (or univariate) probability (or percentile) from the function SDI
#' 
#' 
#' @param (X, Y): The vector of a monthly hydro-climatic variable of n years. ts is the accumulated time scale.
#' @return The multivariate drought index of different time scale from the marginal probability (or percentile) 
#' @export
#' @examples
#X=runif(120, min = 0, max = 100)

#X=c(74.1 ,135.3 , 69.3  ,21.1  ,12.4 , 69.8 , 58.9  ,58.4  ,50.3  ,23.4 , 49.6 ,163.3  ,31.7 , 16.8, 108.7 , 30.8  ,23.7 , 31.1  ,34.2  ,19.0 , 11.9 , 26.5 ,128.8,  79.6  ,30.4 , 35.9 , 45.0  , 7.4  ,23.5 ,109.0 , 23.8  ,87.9 , 22.8  ,54.7  ,39.1 , 80.7  ,36.2 , 92.6,  22.9 ,170.7 , 91.8 ,396.1 ,28.1 ,5.8 ,22.9 ,64.2,70.1 ,85.9)

#Y=c(7.1 ,13.3 , 6.3  ,21.1  ,2.4 , 6.8 , 58.9  ,58.4  ,50.3  ,23.4 , 49.6 ,163.3  ,31.7 , 16.8, 108.7 , 30.8  ,23.7 , 31.1  ,34.2  ,19.0 , 11.9 , 26.5 ,128.8,  79.6  ,30.4 , 35.9 , 45.0  , 7.4  ,23.5 ,109.0 , 23.8  ,87.9 , 22.8  ,54.7  ,39.1 , 80.7  ,36.2 , 92.6,  22.9 ,170.7 , 91.8 ,396.1 ,28.1 ,5.8 ,22.9 ,64.2,70.1 ,85.9)

#' fit<-MDI(X,Y,ts=3) # Compute the 3 month drought index
#' fit$ProbEmp2 #Get the empirival drought index
MDI<-function (X,Y,ts=6) 
{
   
  
  X=matrix(X,ncol=1)
  Y=matrix(Y,ncol=1)
  # Obtain the accumulation of the variable for a specific time scale (ts) 
  
  XA=ACCU(X,ts)
  YA=ACCU(Y,ts)
      
      #  define the SI=Y0 for each month
      #  define the empirical distribution
      
  P_emp=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  P_gam=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)

  Emp2=matrix(NA,nrow=length(X)/12,ncol=12,byrow=T)
  
  
  for  (k in 1:12)
  {
  
    mdx=XA[,k]
    mdy=YA[,k]
    
    # if the first number is NA, exclude the first number and then compute SPI
    
    if (is.na(mdx[1])==TRUE) 
    {
      xd=mdx[2:length(mdx)]
      yd=mdy[2:length(mdy)]
    

      Emp2[2:length(mdx),k]=BiEmp(xd,yd)
      
    }
    
    else
      # if the first number is not NA, take the whole month to compute SPI
    {  
      
      xd=mdx
      yd=mdy    
       
      Emp2[,k]=BiEmp(xd,yd)
      
    }
    
  }
  
  #result<-list(Probemp=P_emp,Probgam=P_gam,ProbEmp2=Emp2)
   
r=qnorm(Emp2)
  result<-list(ProbEmp2=r)
  
  return(result)
}
