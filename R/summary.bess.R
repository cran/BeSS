summary.bess=function(object, ...)
{
  if(object$family=="bess_gaussian") print(cbind(Df=object$s.list,Dev=object$mse))else
    print(cbind(Df=object$s.list,Dev=object$deviance))
}

