bess.lm = function(x, y, beta0, s, max.steps = 20, factor = NULL,
                   weights = rep(1, nrow(x)), normalize = FALSE)
{
  if(missing(beta0)) beta0 = rep(0, ncol(x))
  if(s>length(beta0))
  {stop("s is too large")}
  if(!is.null(factor)){
    if(is.null(colnames(x))) colnames(x) = paste0("X",1:ncol(x),"g")
    if(!is.data.frame(x)) x = as.data.frame(x)
    x[,factor] = apply(x[,factor,drop = FALSE], 2, function(x){
      x = as.factor(x)
    })
    group = rep(1, ncol(x))
    names(group) = colnames(x)
    group[factor] = apply(x[,factor,drop = FALSE], 2,function(x) {length(unique(x))})-1
    Gi = rep(1:ncol(x), times = group)
    beta0 = rep(beta0, times = group)
    x = model.matrix(~., data = x)[,-1]
    fit = gbess.lm(x, y, Gi, beta0, s = s, max.steps = max.steps,
                   weights = weights, normalize = normalize)
    fit$factor = factor
    return(fit)
  }else{
    if(!is.matrix(x)) x = as.matrix(x)
    if(is.null(colnames(x))) colnames(x) = paste0("X",1:ncol(x))
    vn = dimnames(x)[[2]]
    weights = weights/mean(weights)

    fit = bess_lm(x, y, s, max.steps, beta0, weights, normalize)
    beta_out = fit$beta
    names(beta_out) = vn
    coef0 = fit$coef0
    A = which(beta_out!=0)
    xbest = x[,A]
    colnames(xbest) = vn[A]
    bestmodel = lm(y~xbest, weights = weights)

    return(list(family = "bess_gaussian", beta = beta_out, coef0 = fit$coef0,
                nsample = nrow(x), bestmodel = bestmodel,
                mse = fit$mse, nullmse = fit$nullmse,
                AIC = fit$aic, BIC = fit$bic, GIC = fit$gic,
                max.steps = max.steps))
  }
}



