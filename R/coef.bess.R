coef.bess=function(object, sparse=TRUE, ...)
{
   if(sparse==TRUE)
   {
     beta=object$beta
     beta=Matrix(beta,sparse = TRUE)
     return(beta)
   }else return(object$beta)
}

coef.bess.one=function(object,sparse=TRUE, ...)
{
  if(sparse==TRUE)
  {
    beta=object$beta
    beta=matrix(beta,byrow =TRUE)
    beta=Matrix(beta,sparse = TRUE)
    return(beta)
  }else return(object$beta)
}
