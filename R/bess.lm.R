bess.lm=function(x,y,beta0,s,max.steps=20,normalize=FALSE)
{
  if(missing(beta0)) beta0=rep(0,ncol(x))
  if(s>length(beta0))
  {stop("s is too large")}
  n = dim(x)[1]
  m = dim(x)[2]
  im = inactive = seq(m)
  vn = dimnames(x)[[2]]
  one = rep(1,n)
  beta=beta0
  names(beta) = vn
  xs=x
  ys=y
  if(normalize)
  {
    meanx = drop(one %*% x)/n
    x = scale(x, meanx, FALSE)
    mu = mean(y)
    y = drop(y - mu)

    normx = sqrt(drop(one %*% (x^2)))
    nosignal = normx/sqrt(n) < (.Machine$double.eps)
    if (any(nosignal)) {
      ignores = im[nosignal]
      inactive = im[-ignores]
      normx[nosignal] = (.Machine$double.eps) * sqrt(n)
    }    else ignores = NULL
    names(normx) = NULL
    x = sqrt(n)*scale(x, FALSE, normx)

  }

  fit=bess_lm(x,y,s,max.steps,beta0)

  beta=fit$beta
  xbest=xs[,which(beta!=0)]
  bestmodel=lm(ys~xbest)

  lambda=fit$max_T^2/2
  mse=mean((y-x%*%beta)^2)
  nullmse=mean((y^2))
  aic=n*log(mse)+2*s
  bic=n*log(mse)+log(n)*s
  ebic=n*log(mse)+(log(n)+2*log(m))*s
  if(normalize)
  {
    beta=sqrt(n)*beta/normx
    coef0=mu-sum(beta*meanx)
    return(list(family="bess_gaussian",beta=beta,coef0=coef0,nsample=n,bestmodel=bestmodel,
                lambda=lambda,mse=mse,nullmse=nullmse,AIC=aic,BIC=bic,EBIC=ebic))
  }else return(list(family="bess_gaussian",beta=beta,coef0=0,nsample=n,bestmodel=bestmodel,
                    lambda=lambda,mse=mse,nullmse=nullmse,AIC=aic,BIC=bic,EBIC=ebic))
}



