\name{BeSS-internal}
\title{Internal BeSS functions}
\alias{plot_solution}
\alias{plot_loss}
\alias{bess_lm}
\alias{gbess_lm}
\alias{get_A}
\alias{getcox_A}
\alias{gget_A}
\alias{bess.lm}
\alias{bess.glm}
\alias{bess.cox}
\alias{gbess.lm}
\alias{gbess.glm}
\alias{EigenR}
\alias{uniqueR}
\description{Internal BeSS functions}
\usage{
bess.lm(x, y, beta0, s, max.steps=20, factor = NULL,
        weights=rep(1,nrow(x)), normalize=FALSE)
bess.glm(x, y, beta0, intercept=0, s, max.steps=20,
         glm.max=1e6, factor = NULL,
         weights=rep(1,nrow(x)), normalize=FALSE)
bess.cox(x, y, beta0, s, cox.max=20, max.steps=20, weights=rep(1,nrow(x)),
         normalize=FALSE)
gbess.lm(x, y, Gi, beta0, s, max.steps = 20,
         weights=rep(1,nrow(x)), normalize=FALSE)
gbess.glm(x, y, Gi, beta0, intercept=0, s, max.steps = 5, glm.max=1e6,
          weights=rep(1,nrow(x)), normalize=FALSE)		 
}
\author{
  Canhong Wen, Aijun Zhang, Shijie Quan, and Xueqin Wang.
}
\references{
 Wen, C., Zhang, A., Quan, S. and Wang, X. (2017). BeSS: an R package for best subset selection in linear, logistic and CoxPH models. arXiv: 1709.06254.
}
\details{These are not intended for use by users. \code{bess.lm}fit a linear regression model. \code{bess.glm}fit a logistic model. \code{bess.cox}fit a cox model.}
\keyword{internal}
