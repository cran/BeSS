bess.cox = function(x, y, beta0, s, max.steps = 20,
                    cox.max = 20, eta = 0.2, factor = NULL,
                    weights = rep(1, nrow(x)), normalize = FALSE)
{
  if(missing(beta0)) beta0 = rep(0, ncol(x))
  if(is.matrix(y) != 0) y = as.matrix(y)
  if(ncol(y) != 2) stop("Please input y with two columns!")
  if(s > length(beta0))
  {stop("s is too large")}
  if(is.null(colnames(x))) colnames(x) = paste0("X", 1:ncol(x))
  if(!is.null(factor)){
    if(!is.data.frame(x)) x = as.data.frame(x)
    x[,factor] = apply(x[,factor,drop = FALSE], 2, function(x){
      return(as.factor(x))
    })
    group = rep(1, ncol(x))
    names(group) = colnames(x)
    group[factor] = apply(x[,factor,drop = FALSE], 2,function(x) {length(unique(x))})-1
    Gi = rep(1:ncol(x), times = group)
    beta0 = rep(beta0, times = group)
    x = model.matrix(~., data = x)[,-1]
    fit = gbess.cox(x, y, Gi, beta0=beta0, s = s,
                    max.steps = max.steps, cox.max = cox.max,
                    weights = weights, normalize = normalize)
    fit$factor = factor
    return(fit)
  }else{
    x = as.matrix(x)
    n = dim(x)[1]
    p = dim(x)[2]
    vn = dimnames(x)[[2]]
    mark = order(y[,1])
    y = y[mark,]
    x = x[mark,]
    weights = weights[mark]
    weights = weights/mean(weights)
    res = bess_cox(x, y[,2,drop = TRUE], s, max.steps, beta0, weights,
                   cox.max, eta, normalize)
    l = res$l
    A = res$A
    xbest = x[,A]
    #beta = res$beta
    beta = rep(0, p)
    names(beta) = vn
    suppressWarnings(bestmodel <- coxph(Surv(y[,1], y[,2])~xbest, weights = weights, iter.max = cox.max))
    beta[A] = bestmodel$coefficients
    nulldev = -2*bestmodel$loglik[1]

    if(l>max.steps) warning("algorithm did not converge")
    return(list(family = "bess_cox", beta = beta, nsample = n, bestmodel = bestmodel,
                deviance = -2*bestmodel$loglik[2], nulldeviance = nulldev,
                AIC = res$aic, BIC = res$bic, GIC = res$gic,
                max.steps = max.steps))
  }
}
