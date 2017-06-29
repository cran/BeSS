\name{coef.bess}
\alias{coef.bess}
\alias{coef.bess.one}
\title{Provides estimated coefficients from a fitted "bess" object.}
\description{
  Similar to other prediction methods, this function provides estimated coefficients from a fitted "\code{bess}" object or a fitted "\code{bess.one}" object.
}
\usage{
    \method{coef}{bess}(object, sparse=TRUE, type = c("ALL", "AIC", "BIC", "EBIC"),...)
    \method{coef}{bess.one}(object, sparse = TRUE , ...)
}
\arguments{
  \item{object}{A "\code{bess}" project or a "\code{bess.one}" project.}
  \item{sparse}{Logical or NULL, specifying whether the coefficients should be presented as sparse matrix or not.}
  \item{type}{Types of coefficients returned. \code{ type = "AIC"} cooresponds to the coefficient with optimal AIC value; \code{ type = "BIC"} cooresponds to the coefficient with optimal BIC value; \code{ type = "EBIC"} cooresponds to the coefficient with optimal EBIC value; \code{ type = "ALL"} cooresponds to all coefficients in the \code{bess} object. Default is \code{ALL}. }
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
