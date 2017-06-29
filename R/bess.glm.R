bess.glm=function(x,y,beta0,
                  intercept=0,s,
                  max.steps=20,
                  glm.max=1e6,
                  normalize=FALSE)
{
  if(length(unique(y))!=2)  stop("Please input binary variable!")

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



      if(normalize==TRUE)
      {
        meanx = drop(one %*% x)/n
        x = scale(x, meanx, FALSE)

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
      setA=get_A(x,y,beta,intercept,s,rep(0,m))
      p=setA$p
      l=1

      A=list()
      I=list()

      A[[1]]=0
      I[[1]]=seq(m)
      A[[l+1]] = setA$A
      I[[l+1]] = setA$I
      S=1:nrow(x)

      while (l <= max.steps)
      {

        beta[I[[l+1]]] = 0

        if(s>=2)
         {
           logit=glmnet(x[,A[[l+1]]],y,family="binomial",lambda = 0,maxit=glm.max)
           beta[A[[l+1]]]=logit$beta
           coef0=logit$a0
         }else
         {
           logit=glm(y~x[,A[[l+1]]],family="binomial")
           beta[A[[l+1]]]=logit$coefficients[-1]
           coef0=logit$coefficients[1]
         }



        setA=get_A(x,y,beta,coef0,s,A[[l+1]])
        p=setA$p

        A[[l+2]] = setA$A
        I[[l+2]] = setA$I

        if(setequal(A[[l+2]],A[[l]])|setequal(A[[l+2]],A[[l+1]])) {break}
        else{l=l+1
        gc()}
      }
      #dev=logit$deviance
      dev=-2*sum((y*log(p) + (1-y)*log(1-p))[which(p>1e-20&p<1-1e-20)])
      nulldev=-2*sum(y*log(0.5) + (1-y)*log(0.5))
      aic=dev+2*s
      bic=dev+log(n)*s
      ebic=dev+(log(n)+2*log(m))*s

      if(normalize)
      {
        beta=sqrt(n)*beta/normx
        coef0=coef0-sum(beta*meanx)
      }

      return(list(family="bess_binomial",beta=beta,coef0=coef0,
                  deviance=dev,nulldeviance=nulldev,
                  lambda=setA$max_T^2/2,p=p,AIC=aic,BIC=bic,EBIC=ebic))
}





