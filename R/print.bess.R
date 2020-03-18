print.bess=function(x, ...)
{
  if(x$family=="bess_gaussian") print(cbind(Df=x$s.list,Dev=x$mse,AIC=x$AIC,BIC=x$BIC))else
    print(cbind(Df=x$s.list,Dev=x$deviance,AIC=x$AIC,BIC=x$BIC))
}

