summary.bess=function(object, ...)
{
  if(object$family=="bess_gaussian") print(cbind(Df=object$s.list,Dev=object$mse,AIC=object$AIC,BIC=object$BIC))else
    print(cbind(Df=object$s.list,Dev=object$deviance,AIC=object$AIC,BIC=object$BIC))
}

