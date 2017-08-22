\name{print.bess}
\alias{print.bess}
\title{print method for a "bess" object}
\description{
  Print a summary of the "\code{bess}" path at each step along the path.
}
\usage{
  \method{print}{bess}(x, ...)
}
\arguments{
  \item{x}{a "\code{bess}" object}
  \item{\dots}{additional print arguments}
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
print(fit)

}
