predict.bess=function(object, newdata, type = c("ALL", "AIC","BIC"),...)
{
  type <- match.arg(type)
  if(object$family == "bess_gaussian")
  {
    newx = newdata
    betas = object$beta
    coef0 = object$coef0
    y = t(newx%*%betas)+coef0
    if(type == "ALL"){
      return(y)
    }else return(y[which.min(object[[type]]),,drop = TRUE])
  }
  if(object$family == "bess_binomial")
  {
    newx = newdata
    betas = object$beta
    coef = object$coef0
    class = matrix(0,ncol(betas),nrow(newx))
    for(i in 1:ncol(betas))
    {
      class[i,] = as.numeric(exp(newx%*%betas[,i]+coef[i])/(1+exp(newx%*%betas[,i]+coef[i]))>0.5)
      class[i,][is.na(class[i,])] = 1
      if(!is.null(object$y_names))
      {
       class[which(class == 0,arr.ind = T)] = object$y_names[1]
       class[which(class == 1,arr.ind = T)] = object$y_names[2]
      }
    }
    if(type == "ALL"){
      return(class)
    }else return(class[which.min(object[[type]]),,drop = TRUE])
  }
  if(object$family=="bess_cox")
  {
    newx = as.matrix(newdata)
    betas = object$beta

    betax = newx%*%betas
    if(type == "ALL"){
      return(t(betax))
    }else return(betax[,which.min(object[[type]]),drop = TRUE])
  }

}



predict.bess.one=function(object,newdata, ...)
{
  if(object$family=="bess_gaussian")
  {
    newx = newdata
    betas = object$beta
    coef0 = object$coef0

    y = drop(newx %*% betas)+coef0
    return(y)
  }
  if(object$family == "bess_binomial")
  {
    newx = newdata
    betas = object$beta
    coef = object$coef0

    class = as.numeric(exp(newx%*%betas+coef)/(1+exp(newx%*%betas+coef))>0.5)
    class[is.na(class)] = 1
    if(!is.null(object$y_names))
    {
      class[which(class == 0,arr.ind = T)] = object$y_names[1]
      class[which(class == 1,arr.ind = T)] = object$y_names[2]
    }

    return(class)
  }
  if(object$family == "bess_cox")
  {
    newx = as.matrix(newdata)
    betas = object$beta

    betax = newx%*%betas;
    return(drop(betax))
  }

}

