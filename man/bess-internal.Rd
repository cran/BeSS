\name{BeSS-internal}
\title{Internal BeSS functions}
\alias{plot_solution}
\alias{plot_loss}
\alias{abess_lm}
\alias{bess_lm}
\alias{get_A}
\alias{get_A2}
\alias{getcox_A}
\alias{getcox_A2}
\alias{bess.lm}
\alias{bess.glm}
\alias{bess.cox}
\description{Internal BeSS functions}
\usage{
bess.lm(x, y, beta0, s, max.steps=20, normalize=FALSE)
bess.glm(x, y, beta0, intercept=0, s, max.steps=20, glm.max=1e6, normalize=FALSE)
bess.cox(x, y, beta0, s, cox.max=20, max.steps=20, normalize=FALSE)
}
\author{
  Canhong Wen, Aijun Zhang, Shijie Quan, and Xueqin Wang.
}
\details{These are not intended for use by users. \code{bess.lm}fit a linear regression model. \code{bess.glm}fit a logistic model. \code{bess.cox}fit a cox model.}
\keyword{internal}
