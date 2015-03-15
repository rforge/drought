BiEmp=function(X,Y)
  
{
X=matrix(X,ncol=1)
Y=matrix(Y,ncol=1)

n=length(X)

Z=matrix(NA,nrow=length(X),ncol=1)
 
for  (k in 1:n)
{  
  Z[k]=sum((X<=X[k])&(Y<=Y[k])) 
  
  Z[k]=Z[k]/(n+1)
  
} 

return(Z)

}