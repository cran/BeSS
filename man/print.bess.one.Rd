\name{print.bess.one}
\alias{print.bess.one}
\title{print method for a "bess.one" object}
\description{
  Print the primary elements of the "\code{bess.one}" object.
}
\usage{
  \method{print}{bess.one}(x, ...)
}
\arguments{
  \item{x}{a "\code{bess.one}" object}
  \item{\dots}{additional print arguments}
}

\author{
  Canhong Wen, Aijun Zhang, Shijie Quan, and Xueqin Wang.
}
\references{
 Wen, C., Zhang, A., Quan, S. and Wang, X. (2017). BeSS: an R package for best subset selection in linear, logistic and CoxPH models. arXiv: 1709.06254.
}

\seealso{
\code{\link{bess}}, \code{\link{bess.one}}
}

\examples{
data <- gen.data(500, 20, family = "gaussian", 10, 0.2, 1)
fit <- bess.one(data$x, data$y, s = 10, family = "gaussian")
print(fit)

}
