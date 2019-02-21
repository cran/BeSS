summary.bess = function(object, ...){
  beta = object$beta
  if(object$method == "sequential"){
   K.opt.aic = which.min(object$AIC)
   K.opt.bic = which.min(object$BIC)
   K.opt.gic = which.min(object$GIC)
   predictors.aic = beta[, K.opt.aic]
   predictors.bic = beta[, K.opt.bic]
   predictors.gic = beta[, K.opt.gic]
   if(sum(predictors.aic != 0) > 1) predictor.a = "predictors" else predictor.a = "predictor"
   if(sum(predictors.bic != 0) > 1) predictor.b = "predictors" else predictor.b = "predictor"
   if(sum(predictors.gic != 0) > 1) predictor.g = "predictors" else predictor.g = "predictor"
   cat("-------------------------------------------------------------------------------\n")
   cat("    Primal-dual active algorithm with tuning parameter determined by sequential method", "\n\n")
   cat("    Best model determined by AIC includes" , sum(predictors.aic != 0), predictor.a, "with AIC =",
       object$AIC[K.opt.aic], "\n\n")
   cat("    Best model determined by BIC includes" , sum(predictors.bic != 0), predictor.b, "with BIC =",
       object$BIC[K.opt.bic], "\n\n")
   cat("    Best model determined by GIC includes" , sum(predictors.gic != 0), predictor.g, "with GIC =",
       object$GIC[K.opt.gic], "\n")
   cat("-------------------------------------------------------------------------------\n")
  } else {
    cat("------------------------------------------------------------------------------\n")
    cat("    Primal-dual active algorithm with tuning parameter determined by gsection method", "\n\n")
    if(sum(beta[,ncol(beta)] != 0) > 0) predictor = "predictors" else predictor = "predictor"
    cat("    Best model includes", sum(beta[,ncol(beta)] != 0), predictor, "with", "\n\n")
    if(logLik(object)[length(logLik(object))] >= 0)
      cat("    log-likelihood:   ", logLik(object)[length(logLik(object))],"\n") else cat("    log-likelihood:  ", logLik(object)[length(logLik(object))],"\n")

    if(deviance(object)[length(deviance(object))] >= 0)
      cat("    deviance:         ", deviance(object)[length(deviance(object))],"\n") else  cat("    deviance:       ", deviance(object)[length(deviance(object))],"\n")

    if(object$AIC[length(object$AIC)] >= 0)
      cat("    AIC:              ", object$AIC[length(object$AIC)],"\n") else   cat("    AIC:             ", object$AIC[length(object$AIC)],"\n")

    if(object$BIC[length(object$BIC)] >= 0)
      cat("    BIC:              ", object$BIC[length(object$BIC)],"\n") else   cat("    BIC:             ", object$BIC[length(object$BIC)],"\n")

    if(object$GIC[length(object$GIC)] >= 0)
      cat("    GIC:              ", object$GIC[length(object$GIC)],"\n") else    cat("    GIC:            ", object$GIC[length(object$GIC)],"\n")
    cat("------------------------------------------------------------------------------\n")
  }

}
