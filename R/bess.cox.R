bess.cox=function(x,y,beta0,s,
                  cox.max=20,
                  max.steps=20,
                  normalize=FALSE)
{
  if(missing(beta0)) beta0=rep(0,ncol(x))
  if(is.matrix(y)!=0) y=as.matrix(y)
  if(ncol(y)!=2) stop("Please input y with two columns!")
  if(s>length(beta0))
  {stop("s is too large")}

  n = dim(x)[1]
  m = dim(x)[2]
  im = inactive = seq(m)
  vn = dimnames(x)[[2]]
  one = rep(1,n)
  beta=beta0
  names(beta) = vn

  if(normalize==TRUE)
  {
    mark=order(y[,1],decreasing = FALSE)
    y=y[mark,]
    x=x[mark,]

    one = rep(1, n)
    #center
    meanx = drop(one %*% x)/n
    x = scale(x, meanx, FALSE)
    #normalize
    normx = sqrt(drop(one %*% (x^2)))
    nosignal = normx/sqrt(n) < .Machine$double.eps
    if (any(nosignal)) {
      ignores = im[nosignal]
      inactive = im[-ignores]
      normx[nosignal] = (.Machine$double.eps) * sqrt(n)
    }    else ignores = NULL
    names(normx) = NULL
    x = sqrt(n)*scale(x, FALSE, normx)

  }


  ind=which(y[,2]==0)

  setA=getcox_A(x,y,beta,s,rep(0,m),status=ind)
  l=1

  A=list()
  I=list()

  A[[1]]=0
  I[[1]]=seq(m)
  A[[l+1]] = setA$A
  I[[l+1]] = setA$I

  while ((l <= max.steps))
  {

    beta[I[[l+1]]] = 0

    cox=coxph(Surv(y[,1],y[,2])~x[,A[[l+1]]],eps=1e-8,iter.max=cox.max)
    beta[A[[l+1]]]=cox$coefficients

    setA=getcox_A(x,y,beta,s,A[[l+1]],status=ind)

    A[[l+2]] = setA$A
    I[[l+2]] = setA$I

    if(setequal(A[[l+2]],A[[l]])|setequal(A[[l+2]],A[[l+1]])) {break}
    else{l=l+1
    gc()}
  }
  dev=-2*cox$loglik[2]
  nulldev=-2*cox$loglik[1]
  aic=dev+2*s
  bic=dev+log(n)*s
  ebic=dev+(log(n)+2*log(m))*s

  if(normalize==T)
  {
    beta=sqrt(n)*beta/normx
  }

  return(list(family="bess_cox",beta=beta,deviance=dev,
              nulldeviance=nulldev,lambda=setA$max_T^2/2,AIC=aic,BIC=bic,EBIC=ebic))
}




