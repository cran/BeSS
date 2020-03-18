\name{deviance.bess}
\alias{deviance.bess}
\alias{deviance.bess.one}
\title{Extract the deviance from a "bess" object.}
\description{
  Similar to other deviance methods, which returns deviance from a fitted "\code{bess}" object or a fitted "\code{bess.one}" object.
}
\usage{
  \method{deviance}{bess}(object,...)
  \method{deviance}{bess.one}(object,...)
}
\arguments{
  \item{object}{Output from the \code{bess} function or the \code{bess.one} function.}
  \item{\dots}{Additional arguments affecting the predictions produced.}
}
\value{The value of the deviance extracted from the "bess" object.}

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
deviance(fit)

}
