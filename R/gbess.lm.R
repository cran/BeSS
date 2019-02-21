gbess.lm = function(x, y, Gi, beta0, s, max.steps = 20,
                    weights = rep(1, nrow(x)), normalize = FALSE)
{
  if(missing(beta0)) beta0 = rep(0, ncol(x))
  if(s>length(beta0))
  {stop("s is too large")}
  # initial
  n = dim(x)[1]
  p = dim(x)[2]
  vn = dimnames(x)[[2]]
  beta = beta0
  names(beta) = vn

  weights = weights/mean(weights)
  orderGi = order(Gi)
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
  fit = gbess_lm(X = x, y = y, G = Gi, index = gi_index, orderGi = orderGi, PhiG = PhiG, invPhiG = invPhiG,
                 T0 = s, max_steps = max.steps, beta0 = beta0, weights = weights,
                 n = n, p = p, N = N, normal = normalize)
  beta_tmp = fit$beta
  beta = rep(0, p)
  beta[orderGi] = beta_tmp
  names(beta) = vn
  xbest = x[, which(beta != 0)]
  bestmodel = lm(y~xbest, weights = weights)
  s = length(which(beta != 0))
  mse = mean(weights*(y-x%*%beta)^2)
  nullmse = mean(weights*(y^2))
  aic = n*log(mse)+2*s
  bic = n*log(mse)+log(n)*s
  gic = n*log(mse)+log(p)*log(log(n))*s
  return(list(family = "bess_gaussian", beta = beta, coef0 = fit$coef0, nsample = n, bestmodel = bestmodel,
              mse = mse, nullmse = nullmse, AIC = aic, BIC = bic, GIC = gic, max.steps = max.steps,
              gr_size = fit$gr_size))
}


