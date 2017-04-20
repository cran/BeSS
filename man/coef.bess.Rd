\name{coef.bess}
\alias{coef.bess}
\alias{coef.bess.one}
\title{Provides estimated coefficients from a fitted "bess" object.}
\description{
  Similar to other prediction methods, this function provides estimated coefficients from a fitted "\code{bess}" object or a fitted "\code{bess.one}" object.
}
\usage{
    \method{coef}{bess}(object, sparse = TRUE , ...)
    \method{coef}{bess.one}(object, sparse = TRUE , ...)
}
\arguments{
  \item{object}{A "\code{bess}" project or a "\code{bess.one}" project.}
  \item{sparse}{Logical or NULL, specifying whether the coefficients should be presented as sparse matrix or not.}
  \item{\dots}{Other arguments.}
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
coef(fit, sparse=TRUE)  # The estimated coefficients
}
