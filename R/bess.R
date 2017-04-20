bess = function(x, y, family = c("gaussian", "binomial", "cox"),
                method = "gsection", s.min = 1,
                s.max,
                s.list,
                K.max = 20,
                max.steps = 15,
                glm.max = 1e6,
                cox.max = 20,
                epsilon = 1e-4,
                normalize = FALSE)
{
  family <- match.arg(family)
  if(ncol(x)==1|is.vector(x)) stop("x should be two columns at least!")
  if(missing(family)) stop("Please input family!")
  if(family=="binomial")
  {
    if(length(unique(y))!=2)  stop("Please input binary variable!")else
      if(setequal(y_names<-unique(y),c(0,1))==FALSE)
      {
        y[which(y==unique(y)[1])]=0
        y[which(y==unique(y)[2])]=1
        y=as.numeric(y)
      }
  }
  if(family=="cox")
  {
    if(!is.matrix(y)) y=as.matrix(y)
    if(ncol(y)!=2) stop("Please input y with two columns!")
  }
  if(is.vector(y))
  {
    if(nrow(x)!=length(y)) stop("Rows of x must be the same as length of y!")
  }else{
    if(nrow(x)!=nrow(y)) stop("Rows of x must be the same as rows of y!")
  }
  if(method=="sequential"&missing(s.list)) s.list=1:min(ncol(x),round(nrow(x)/log(nrow(x))))


  if(missing(s.max)) s.max=min(ncol(x),round(nrow(x)/log(nrow(x))))

#normalize x
  x=as.matrix(x)
  nm = dim(x)
  n = nm[1]
  m = nm[2]
  im = inactive = seq(m)
  one = rep(1, n)
  vn = dimnames(x)[[2]]

  meanx = drop(one %*% x)/n
  x = scale(x, meanx, FALSE)



  if(family=="gaussian")
  {
    mu = mean(y)
    y = drop(y - mu)

    normx = sqrt(drop(one %*% (x^2)))
    nosignal = normx/sqrt(n) < .Machine$double.eps
    if (any(nosignal)) {
      ignores = im[nosignal]
      inactive = im[-ignores]
      normx[nosignal] = (.Machine$double.eps) * sqrt(n)
    }    else ignores = NULL
    names(normx) = NULL
    x = sqrt(n)*scale(x, FALSE, normx)

    #initial beta
    beta0=rep(0,m)
    gc()

    if(method=="gsection")
    {
     fit=abess_lm(x,y,s.min,s.max,K.max,max.steps,beta0,epsilon)
     beta.fit=fit$beta[,1:fit$k]
     mse=fit$mse[1:fit$k]
     lambda=fit$lambda[1:fit$k]^2/2
     nullmse=fit$nullmse
     s.list=fit$s.list[1:fit$k]
    }
    if(method=="sequential")
    {
      fit_L=bess.lm(x,y,s=s.list[length(s.list)],max.steps = max.steps,beta0=beta0)
      mse_L=fit_L$mse
      #cat(mse_L,"\n")
      nullmse=sum(y^2)/(2*n)
      #cat(nullmse,"\n")
      beta.fit=matrix(beta0,m,1)
      mse=vector()
      lambda=vector()
      for(k in 1:length(s.list))
      {
        cat("select",s.list[k],"variables","\n")
        fit=bess.lm(x,y,s=s.list[k],max.steps = max.steps,beta0=beta.fit[,k,drop=TRUE])
        #estimator of beta for T=k*Tao
        beta.fit=cbind(beta.fit,fit$beta)
        #mse
        mse[k]=fit$mse
        #lambda
        lambda[k]=fit$lambda

        #cat(abs((mse[k]-mse_L)/((nullmse-mse_L)*(L-s.list[k]))),"\n")
        if(s.list[k]==s.list[length(s.list)]) break
        if(abs((mse[k]-mse_L)/((nullmse-mse_L)*(s.list[length(s.list)]-s.list[k])))<epsilon) break
      }
      beta.fit=beta.fit[,-1,drop=FALSE]
      s.list=s.list[1:k]
    }


    beta.fit=sqrt(n)*(beta.fit)/normx
    coef0=mu-sum(beta.fit*meanx)

    out=list(family="bess_gaussian",beta=beta.fit,coef0=coef0,
             s.list=s.list,meanx=meanx,normx=normx,meany=mu,
             mse=mse,nullmse=nullmse,lambda=lambda)
    class(out)="bess"
    return(out)
  }

  if(family=="binomial")
  {
    normx = sqrt(drop(one %*% (x^2)))
    nosignal = normx/sqrt(n) < .Machine$double.eps
    if (any(nosignal)) {
      ignores = im[nosignal]
      inactive = im[-ignores]
      normx[nosignal] = (.Machine$double.eps) * sqrt(n)
    }    else ignores = NULL
    names(normx) = NULL
    x = sqrt(n)*scale(x, FALSE, normx)

    beta0=rep(0,m)
    intercept=0
    gc()

    if(method=="gsection")
    {
    k = 1
    sL=s.min
    sR=s.max
    beta0R=beta0
    beta0L=beta0
    coef0L=intercept
    coef0R=intercept

    fit_L1=bess.glm(x=x,y=y,
                    beta0=beta0L,
                    intercept=coef0L,
                    s=s.list[length(s.list)],
                    glm.max=glm.max,
                    max.steps=max.steps)

    nulldev=fit_L1$nulldeviance

    fit_L=fit_L1
    fit_R=bess.glm( x=x,y=y,
                    beta0=beta0R,
                    intercept=coef0R,
                    s=sR,
                    glm.max=glm.max,
                    max.steps=max.steps)

    beta.fit=cbind(fit_L$beta,fit_R$beta)
    coef0.fit=c(fit_L$coef0,fit_R$coef0)
    dev=c(fit_L$deviance,fit_R$deviance)
    lambda=c(fit_L$lambda,fit_R$lambda)

    beta0M=fit_R$beta
    coef0M=fit_R$coef0
    s.list=c(sL,sR)
    while(k<=K.max)
    {
      sM <- round(sL + (sR-sL)*0.618)
      s.list=c(s.list,sM)
      fit_M=bess.glm(x=x,y=y,
                     beta0=beta0M,
                     intercept=coef0M,
                     s=sM,
                     glm.max=glm.max,
                     max.steps=max.steps)
      cat(k,"-th iteration s.left:",sL," s.split:",sM," s.right:",sR,"\n",sep="")

      beta0M=fit_M$beta
      beta.fit=cbind(beta.fit,beta0M)

      coef0M=fit_M$coef0
      coef0.fit=c(coef0.fit,coef0M)

      dev=c(dev,fit_M$deviance)
      lambda=c(lambda,fit_M$lambda)
      if(abs(fit_L$deviance-fit_M$deviance)/abs(nulldev*(sM-sL)) > epsilon &
         abs(fit_R$deviance-fit_M$deviance)/abs(nulldev*(sM-sR)) < epsilon)
      {
        sR <- sM
        fit_R=fit_M
      }else if(abs(fit_L$deviance-fit_M$deviance)/abs(nulldev*(sM-sL)) > epsilon &
               abs(fit_R$deviance-fit_M$deviance)/abs(nulldev*(sM-sR)) > epsilon)
      {
        sL <- sM
        fit_L=fit_M
      }else
      {
        sR=sM
        fit_R=fit_M
        sL=s.min
        fit_L=fit_L1
      }

      if(sR-sL==1) break
      fit_ML=bess.glm(x=x,y=y,
                     beta0=beta0M,
                     intercept=coef0M,
                     s=sM,
                     glm.max=glm.max,
                     max.steps=max.steps)

      fit_MR=bess.glm(x=x,y=y,
                     beta0=beta0M,
                     intercept=coef0M,
                     s=sM,
                     glm.max=glm.max,
                     max.steps=max.steps)
      #if(abs(fit_ML$deviance-fit_M$deviance)/abs(fit_M$deviance) > epsilon &
        # abs(fit_MR$deviance-fit_M$deviance)/abs(fit_M$deviance) < epsilon)
     # {break}
      if(abs(fit_ML$deviance-fit_M$deviance)/abs(nulldev) > epsilon &
         abs(fit_MR$deviance-fit_M$deviance)/abs(nulldev) < epsilon)
      {break}

    k=k+1
    }
    }

    if(method=="sequential")
    {
     fit_L=bess.glm(x=x,y=y,beta0=beta0,intercept=intercept,
                     s=s.list[length(s.list)],glm.max=glm.max,max.steps = max.steps)
      beta.fit=matrix(beta0,m,1)
      coef0.fit=intercept
      dev_L=fit_L$deviance
      nulldev=-2*sum((y*log(0.5) + (1-y)*log(0.5)))

      if(abs(dev_L/nulldev)>0.5) dev_L=0
      dev=vector()
      lambda=vector()
      for(k in 1:length(s.list))
      {
        cat("select",s.list[k],"variables","\n")

        fit=bess.glm(x=x,y=y,beta0=beta.fit[,k],intercept=coef0.fit[k],
                          s=s.list[k],glm.max=glm.max,
                          max.steps=max.steps)

        dev[k]=fit$deviance

        lambda[k]=fit$lambda
        beta.fit=cbind(beta.fit,fit$beta)
        coef0.fit=c(coef0.fit,fit$coef0)
        #cat(abs((dev[k]-dev_L)/((s.list[length(s.list)]-s.list[k])*(nulldev))),"\n")
        if(s.list[k]==s.list[length(s.list)]) break
        if(abs((dev[k]-dev_L)/((s.list[length(s.list)]-s.list[k])*(nulldev-dev_L)))<epsilon) break

    }
      beta.fit=beta.fit[,-1,drop=FALSE]
      coef0.fit=coef0.fit[-1]
      s.list=s.list[1:k]
    }

    beta.fit=sqrt(n)*(beta.fit)/normx
    coef0.fit=coef0.fit-drop(t(beta.fit)%*%meanx)


    if(!setequal(y_names,c(0,1)))
    {
      out=list(family="bess_binomial",beta=beta.fit,coef0=coef0.fit,s.list=s.list,
         meanx=meanx,normx=normx,dev=dev,nulldev=nulldev,lambda=lambda,y_names=y_names)
      class(out)="bess"
      return(out)
    }else
    {
      out=list(family="bess_binomial",beta=beta.fit,coef0=coef0.fit,s.list=s.list,
               meanx=meanx,normx=normx,deviance=dev,nulldeviance=nulldev,lambda=lambda)
      class(out)="bess"
      return(out)
    }
  }

  if(family=="cox")
  {
    #normalize
    normx = sqrt(drop(one %*% (x^2)))
    mark=order(y[,1],decreasing = FALSE)
    y=y[mark,]

    nosignal = normx/sqrt(n) < .Machine$double.eps
    if (any(nosignal)) {
      ignores = im[nosignal]
      inactive = im[-ignores]
      normx[nosignal] = (.Machine$double.eps) * sqrt(n)
    }    else ignores = NULL
    names(normx) = NULL
    x = sqrt(n)*scale(x, FALSE, normx)
    x=x[mark,]

    beta0=rep(0,m)
    gc()

    if(method=="gsection")
    {
    k = 1
    sL=s.min
    sR=s.max
    beta0R=beta0
    beta0L=beta0

    fit_L1=bess.cox(x,y,
                    beta0=beta0L,
                    s=sL,
                    cox.max=cox.max,
                    max.steps=max.steps)

    nulldev=fit_L1$nulldeviance
    fit_L=fit_L1
    fit_R=bess.cox(x,y,
                   beta0=beta0R,
                   s=sR,
                   cox.max=cox.max,
                   max.steps=max.steps)

    beta.fit=cbind(fit_L$beta,fit_R$beta)
    dev=c(fit_L$deviance,fit_R$deviance)
    lambda=c(fit_L$lambda,fit_R$lambda)

    beta0M=fit_R$beta
    s.list=c(sL,sR)

    while(k<=K.max)
    {
      sM <- round(sL + (sR-sL)*0.618)
      s.list=c(s.list,sM)
      fit_M=bess.cox(x,y,
                     beta0=beta0M,
                     s=sM,
                     cox.max=cox.max,
                     max.steps=max.steps)
      cat(k,"-th iteration s.left:",sL," s.split:",sM," s.right:",sR,"\n",sep="")

      beta0M=fit_M$beta
      beta.fit=cbind(beta.fit,beta0M)

      dev=c(dev,fit_M$deviance)
      lambda=c(lambda,fit_M$lambda)
      if(abs(fit_L$deviance-fit_M$deviance)/abs(nulldev*(sM-sL)) > epsilon &
         abs(fit_R$deviance-fit_M$deviance)/abs(nulldev*(sM-sR)) < epsilon)
      {
        sR <- sM
        fit_R=fit_M
      }else if(abs(fit_L$deviance-fit_M$deviance)/abs(nulldev) > epsilon &
               abs(fit_R$deviance-fit_M$deviance)/abs(nulldev) > epsilon)
      {
        sL <- sM
        fit_L=fit_M
      }else
      {
        sR=sM
        fit_R=fit_M
        sL=s.min
        fit_L=fit_L1
      }

      if(sR-sL==1) break
      fit_ML=bess.cox(x,y,
                      beta0=beta0M,
                      s=sM,
                      cox.max=cox.max,
                      max.steps=max.steps)

      fit_MR=bess.cox(x,y,
                      beta0=beta0M,
                      s=sM,
                      cox.max=cox.max,
                      max.steps=max.steps)
      if(abs(fit_ML$deviance-fit_M$deviance)/abs(fit_M$deviance) > epsilon &
         abs(fit_MR$deviance-fit_M$deviance)/abs(fit_M$deviance) < epsilon)
      {break}

      k=k+1
    }
  }

    if(method=="sequential")
    {
      fit_L=bess.cox(x=x,y=y,beta0=beta0,s=s.list[length(s.list)],cox.max=cox.max,
                     max.steps=max.steps)
      beta.fit=matrix(beta0,m,1)
      dev_L=fit_L$deviance
     # cat(dev_L,"\n")
      nulldev=fit_L$nulldeviance

      dev=vector()
      lambda=vector()

      for(k in 1:length(s.list))
      {
        cat("select",s.list[k],"variables","\n")

        fit=bess.cox(x=x,y=y,beta0=beta.fit[,k],
                          s=s.list[k],cox.max=cox.max,
                          max.steps=max.steps)

        dev[k]=fit$deviance
        #cat(dev[k],"\n")
        lambda[k]=fit$lambda
        beta.fit=cbind(beta.fit,fit$beta)


        #cat(abs((dev[k]-dev_L)/((nulldev-dev_L)*(s.list[length(s.list)]-s.list[k]))),"\n")
        if(s.list[k]==s.list[length(s.list)]) break
        if(abs((dev[k]-dev_L)/((nulldev-dev_L)*(s.list[length(s.list)]-s.list[k])))<epsilon) break
      }

      beta.fit=beta.fit[,-1,drop=FALSE]
      s.list=s.list[1:k]
    }

    beta.fit=sqrt(n)*(beta.fit)/normx

    out=list(family="bess_cox",beta=beta.fit,s.list=s.list,meanx=meanx,
             normx=normx,deviance=dev,nulldeviance=nulldev,lambda=lambda)
    class(out)="bess"
    return(out)
  }

}


