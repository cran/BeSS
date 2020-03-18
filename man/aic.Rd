\name{aic}
\alias{aic}
\alias{bic}
\alias{ebic}
\title{Extract the IC from a "bess" object.}
\description{
  These functions are used by bess to compute Information Criteria from a fitted model object.
}
\usage{
  aic(object,...)
  bic(object,...)
  ebic(object,...)
}
\arguments{
  \item{object}{Output from the \code{bess} function or the \code{bess.one} function.}
  \item{\dots}{Additional arguments affecting the predictions produced.}
}
\value{The value of Information Criteria extracted from the "bess" object.}

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
aic(fit)
bic(fit)
ebic(fit)

}
