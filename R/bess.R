bess = function(x, y, family = c("gaussian", "binomial", "cox"),
                method = "gsection", s.min = 1,
                s.max,
                s.list,
                K.max = 20,
                max.steps = 20,
                glm.max = 20,
                cox.max = 20,
                eta = 0.2,
                factor = NULL,
                epsilon = 1e-4,
                weights = rep(1,nrow(x)),
                ic.type = c("GIC", "AIC", "BIC"),
                warm.start = TRUE,
                normalize = TRUE)
{
  family = match.arg(family)
  ic.type = match.arg(ic.type)
  if(ncol(x) == 1|is.vector(x)) stop("x should be two columns at least!")
  if(method == "sequential"&missing(s.list)) s.list = 1:min(ncol(x), round(nrow(x)/log(nrow(x))))
  if(missing(s.max)) s.max = min(ncol(x),round(nrow(x)/log(nrow(x))))
  if(missing(family)) stop("Please input family!")
  if(!is.null(factor)) method = "sequential"
  if(family == "binomial")
  {
    if(is.factor(y)){
      y = as.character(y)
    }
    if(length(unique(y)) != 2)  stop("Please input binary variable!")else
      if(setequal(y_names<-unique(y), c(0,1)) == FALSE)
      {
        y[which(y == unique(y)[1])] = 0
        y[which(y == unique(y)[2])] = 1
        y = as.numeric(y)
      }
  }
  if(family == "cox")
  {
    if(!is.matrix(y)) y = as.matrix(y)
    if(ncol(y) != 2) stop("Please input y with two columns!")
  }
  if(is.vector(y))
  {
    if(nrow(x) != length(y)) stop("Rows of x must be the same as length of y!")
  }else{
    if(nrow(x) != nrow(y)) stop("Rows of x must be the same as rows of y!")
  }
  if(!is.null(factor)){
    if(is.null(colnames(x))) colnames(x) = paste0("X",1:ncol(x),"g")
  }else{
    if(is.null(colnames(x))) colnames(x) = paste0("X",1:ncol(x))
  }

  if(!is.matrix(x)) x = as.matrix(x)
  n = dim(x)[1]
  p = dim(x)[2]
  weights = weights/mean(weights)
  beta0 = rep(0,ncol(x))
  if(!is.null(factor)){
    if(!is.data.frame(x)) x = as.data.frame(x)
    x[,factor] = apply(x[,factor,drop=FALSE], 2, function(x){
      x = as.factor(x)
    })
    group = rep(1, ncol(x))
    names(group) = colnames(x)
    group[factor] = apply(x[,factor,drop=FALSE], 2, function(x) {length(unique(x))})-1
    beta0 = rep(beta0, times = group)
    Gi = rep(1:ncol(x), times = group)
    orderGi = order(Gi)
    x = model.matrix(~., data = x)[,-1]
    p = dim(x)[2]
  }
  vn = colnames(x)
  if(family == "gaussian")
  {
    if(method == "gsection")
    {
      fit = bess_lm_gs(x, y, s.min, s.max, K.max, max.steps,
                       epsilon, beta0, weights, warm.start, normalize)
      beta.fit = fit$beta
      coef0 = fit$coef0
      nullmse = fit$nullmse
      mse = fit$mse
      aic = fit$aic
      bic = fit$bic
      gic = fit$gic
      s.list = fit$s_list
      xbest = x[, which(beta.fit[,ncol(beta.fit), drop=TRUE]!=0)]
      bestmodel = lm(y~xbest, weights = weights)
    }
    if(method == "sequential")
    {
      if(!is.null(factor)){
        x = x[,orderGi]
        Gi = Gi[orderGi]
        gi = unique(Gi)
        gi_index = match(gi, Gi)
        N = length(gi)
        PhiG = lapply(1:N, function(i){
          idx = which(Gi==i)
          if(length(idx) == 1)
            return(-sqrt(t(x[,idx])%*%x[,idx])) else{
              return(-EigenR(t(x[,idx])%*%x[,idx]))
            }
        })
        invPhiG = lapply(PhiG, solve)
        fit = gbess_lms(X = x, y = y, G = Gi, index = gi_index, orderGi = orderGi, PhiG = PhiG, invPhiG = invPhiG,
                        T_list = s.list, max_steps = max.steps, beta0 = beta0, weights = weights,
                        n = n, p = p, N = N, warm_start = warm.start, normal = normalize)
      }else{
        fit = bess_lms(x, y, s.list, max.steps, beta0, weights, warm_start = warm.start, normal = normalize)
      }
      beta.fit = fit$beta
      coef0 = fit$coef0
      nullmse = fit$nullmse
      mse = fit$mse
      aic = fit$aic
      bic = fit$bic
      gic = fit$gic
      if(ic.type == "AIC") s_min = which.min(aic)
      if(ic.type == "BIC") s_min = which.min(bic)
      if(ic.type == "GIC") s_min = which.min(gic)
      xbest = x[,which(beta.fit[,s_min,drop = TRUE]!=0)]
      bestmodel = lm(y~xbest, weights = weights)
    }
    colnames(beta.fit) = s.list
    rownames(beta.fit) = vn
    names(coef0) = s.list
    out = list(family = "bess_gaussian", method = method, beta = beta.fit, coef0 = coef0,
             s.list = s.list, nsample = nrow(x), bestmodel = bestmodel,
             mse = mse, nullmse = nullmse, AIC = aic, BIC = bic, GIC = gic,
             max.steps = max.steps, factor = factor)
    class(out) = "bess"
    return(out)
  }
  if(family == "binomial")
  {
    if(method=="gsection")
    {
      intercept = 0
      fit = bess_glm_gs(x, y, s.min, s.max, K.max, max.steps,
                        epsilon, beta0, intercept, weights,
                        warm.start, normalize, glm.max)
      beta.fit = fit$beta
      coef0.fit = fit$coef0
      nulldev = fit$nulldev
      dev = fit$dev
      aic = fit$aic
      bic = fit$bic
      gic = fit$gic
      s.list = fit$s_list
      xbest = x[,which(beta.fit[,ncol(beta.fit),drop = TRUE]!=0)]
      suppressWarnings(bestmodel <- glm(y~xbest, family=binomial, weights=weights))
    }
    if(method == "sequential")
    {
      if(is.null(factor)) {
        fit = bess_glms(x, y, s.list, max.steps, beta0, 0, weights,
                        warm_start = warm.start, glm_max = glm.max,
                        normal = normalize)

        beta.fit = fit$beta
        coef0.fit = fit$coef0
        dev = fit$deviance
        aic = fit$aic
        bic = fit$bic
        gic = fit$gic
      } else {
        xs = x
        tmp = Normalize2(x, weights)
        x = tmp$X
        meanx = tmp$meanx
        normx = tmp$normx
        rm(tmp)
        intercept = 0
        nulldev = -2*sum(weights*(y*log(0.5) + (1-y)*log(0.5)))
        dev = vector()
        aic = vector()
        bic = vector()
        gic = vector()
        for(k in 1:length(s.list)) {
          fit = gbess.glm(x = x, y = y, Gi = Gi, beta0 = beta0, intercept = intercept,
                          s = s.list[k], max.steps = max.steps,
                          weights = weights)
          if(warm.start) {
            beta0 = fit$beta
            intercept = fit$coef0
          }
          if(k == 1) {
            beta.fit = matrix(fit$beta)
            coef0.fit = fit$coef0
          } else {
            beta.fit = cbind(beta.fit, fit$beta)
            coef0.fit = c(coef0.fit, fit$coef0)
          }
          s.list[k] = fit$gr_size
          dev[k] = fit$deviance
          aic[k] = fit$AIC
          bic[k] = fit$BIC
          gic[k] = fit$GIC
        }
        beta.fit = sqrt(n)*(beta.fit)/normx
        coef0.fit = coef0.fit - drop(t(beta.fit)%*%meanx)
      }
      if(ic.type == "AIC") s_min = which.min(aic)
      if(ic.type == "BIC") s_min = which.min(bic)
      if(ic.type == "GIC") s_min = which.min(gic)
      if(is.null(factor)) {
        xbest = x[,which(beta.fit[,s_min,drop = TRUE]!=0)]
      } else {
        xbest = xs[,which(beta.fit[,s_min,drop = TRUE]!=0)]
      }
    }
    suppressWarnings(bestmodel <- glm(y~xbest, family = binomial, weights = weights))
    colnames(beta.fit) = s.list
    rownames(beta.fit) = vn
    names(coef0.fit) = s.list
    nulldev = -2*sum(weights*(y*log(0.5) + (1-y)*log(0.5)))
    if(!setequal(y_names, c(0,1)))
    {
      out = list(family = "bess_binomial", method = method, beta = beta.fit, coef0 = coef0.fit,
                 s.list = s.list, nsample = nrow(x), bestmodel = bestmodel,
                 deviance = dev, nulldeviance = nulldev, AIC = aic, BIC = bic, GIC = gic,
                 y_names = y_names, max.steps = max.steps, factor = factor)
      class(out) = "bess"
      return(out)
    }else
    {
      out = list(family = "bess_binomial", method = method, beta = beta.fit, coef0 = coef0.fit,
                 s.list = s.list, nsample = nrow(x), bestmodel = bestmodel,
                 deviance = dev, nulldeviance = nulldev, AIC = aic, BIC = bic, GIC = gic,
                 max.steps = max.steps, factor = factor)
      class(out) = "bess"
      return(out)
    }
  }
  if(family == "cox")
  {
    mark = order(y[,1])
    y = y[mark,]
    weights = weights[mark]
    x = x[mark,]
    beta0 = rep(0, p)
    if(method=="gsection")
    {
      fit = bess_cox_gs(x, y[,2,drop = TRUE], s.min, s.max, K.max, max.steps,
                        epsilon, beta0, weights,
                        warm.start, cox.max, eta, normalize)
      beta.fit = fit$beta
      nulldev = fit$nulldev
      dev = fit$dev
      aic = fit$aic
      bic = fit$bic
      gic = fit$gic
      s.list = fit$s_list
      xbest = x[,which(beta.fit[,ncol(beta.fit),drop = TRUE]!=0)]
      suppressWarnings(bestmodel <- coxph(Surv(y[,1],y[,2])~xbest, iter.max = cox.max, weights = weights))
      beta.fit[which(beta.fit[,ncol(beta.fit)]!=0),ncol(beta.fit)] = bestmodel$coefficients
      dev[length(dev)] = bestmodel$loglik[2]*(-2)
    }
    if(method == "sequential")
    {
      if(is.null(factor)){
        fit = bess_coxs(x, y[,2,drop = TRUE], s.list, max.steps, beta0, weights,
                        eta = eta, normal = normalize)
        beta.fit = fit$beta
        dev = fit$deviance
        aic = fit$aic
        bic = fit$bic
        gic = fit$gic
      }else{
        ys = y
        xs = x
        tmp = Normalize2(x, weights)
        x = tmp$X
        meanx = tmp$meanx
        normx = tmp$normx
        rm(tmp)
        dev = vector()
        aic = vector()
        bic = vector()
        gic = vector()
        for(k in 1:length(s.list))
        {
          if(k == 1){
            fit = gbess.cox(x = x, y = y, Gi = Gi, beta0 = beta0,
                            s = s.list[k], cox.max = cox.max,
                            max.steps = max.steps, weights = weights)
            beta.fit = matrix(fit$beta)
            s.list[k] = fit$gr_size
          }else{
            fit = gbess.cox(x, y, Gi = Gi, beta0 = beta.fit[,k-1,drop = TRUE],
                            s = s.list[k], cox.max = cox.max,
                            max.steps = max.steps, weights = weights)
            beta.fit = cbind(beta.fit, fit$beta)
            s.list[k] = fit$gr_size
          }
          dev[k] = fit$deviance
          aic[k] = fit$AIC
          bic[k] = fit$BIC
          gic[k] = fit$GIC
        }
        beta.fit = sqrt(n)*(beta.fit)/normx
      }
      if(ic.type == "AIC") s_min = which.min(aic)
      if(ic.type == "BIC") s_min = which.min(bic)
      if(ic.type == "GIC") s_min = which.min(gic)
      if(is.null(factor)) {
        xbest = x[,which(beta.fit[,s_min,drop = TRUE]!=0)]
        suppressWarnings(bestmodel <- coxph(Surv(y[,1],y[,2])~xbest, iter.max = cox.max, weights = weights))
        beta.fit[which(beta.fit[,s_min]!=0),s_min] = bestmodel$coefficients
      } else {
        xbest = xs[,which(beta.fit[,s_min,drop = TRUE]!=0)]
        suppressWarnings(bestmodel <- coxph(Surv(y[,1],ys[,2])~xbest, iter.max = cox.max, weights = weights))
      }
    }
    colnames(beta.fit) = s.list
    rownames(beta.fit) = vn
    nulldev = bestmodel$loglik[1]*(-2)
    out = list(family = "bess_cox", method = method, beta = beta.fit, s.list = s.list,
               nsample = nrow(x), bestmodel = bestmodel,
               deviance = dev, nulldeviance = nulldev, AIC = aic, BIC = bic, GIC = gic,
               max.steps = max.steps, factor = factor)
    class(out) = "bess"
    return(out)
  }
}
