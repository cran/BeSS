\name{summary.bess}
\alias{summary.bess}
\title{summary method for a "bess" object}
\description{
  Print a summary of the "\code{bess}" path at each step along the path.
}
\usage{
  \method{summary}{bess}(object, ...)
}
\arguments{
  \item{object}{a "\code{bess}" project}
  \item{\dots}{additional summary arguments}
}

\author{
  Canhong Wen, Aijun Zhang, Shijie Quan, and Xueqin Wang.
}
\references{
 Wen, C., Zhang, A., Quan, S., and Wang, X. (2017) BeSS: A R package for best subset selection in GLM and CoxPH Models, Technical reports.
}

\seealso{
\code{\link{bess}}, \code{\link{bess.one}}
}

\examples{
data <- gen.data(500, 20, family = "gaussian", 10, 0.2, 1)
fit <- bess(data$x, data$y, family = "gaussian")
summary.bess(fit)

}
