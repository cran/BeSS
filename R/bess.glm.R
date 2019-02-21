bess.glm = function(x, y, beta0, intercept = 0, s, max.steps = 20,
                    glm.max = 20, factor = NULL,
                    weights = rep(1, nrow(x)), normalize = FALSE)
{
  if(length(unique(y)) != 2)  stop("Please input binary variable!")
  if(missing(beta0)) beta0 = rep(0, ncol(x))
  if(s > length(beta0))
  {stop("s is too large")}
  if(is.null(colnames(x))) colnames(x) = paste0("X",1:ncol(x),"g")
  if(!is.null(factor)){
    if(!is.data.frame(x)) x = as.data.frame(x)
    x[,factor] = apply(x[,factor,drop=FALSE], 2, function(x){
      return(as.factor(x))
    })
    group = rep(1, ncol(x))
    names(group) = colnames(x)
    group[factor] = apply(x[,factor,drop=FALSE], 2,function(x) {length(unique(x))})-1
    Gi = rep(1:ncol(x), times = group)
    beta0 = rep(beta0, times = group)
    x = model.matrix(~., data = x)[,-1]
    fit = gbess.glm(x, y, Gi, beta0=beta0, intercept=intercept, s = s,
                    max.steps = max.steps,
                    weights = weights, normalize = normalize)
    fit$factor = factor
    return(fit)
  }else{
    x = as.matrix(x)
    n = dim(x)[1]
    p = dim(x)[2]
    vn = dimnames(x)[[2]]
    weights = weights/mean(weights)
    res = bess_glm(x, y, s, max.steps, beta0, intercept, weights,
                   glm.max, normalize)
    l = res$l
    A = res$A
    xbest = x[,A]
    beta = res$beta
    names(beta) = vn
    nulldev = -2*sum(weights*(y*log(0.5) + (1-y)*log(0.5)))
    suppressWarnings(bestmodel <- glm(y~xbest, family = binomial, weights = weights))
    if(l>max.steps) warning("algorithm did not converge")
    return(list(family = "bess_binomial", beta = beta, coef0 = res$coef0,
                nsample = n, bestmodel = bestmodel,
                deviance = res$deviance, nulldeviance = nulldev,
                AIC = res$aic, BIC = res$bic, GIC = res$gic,
                max.steps = max.steps))
  }
}





