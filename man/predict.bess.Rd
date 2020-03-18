\name{predict.bess}
\alias{predict.bess}
\alias{predict.bess.one}
\title{make predictions from a "bess" object.}
\description{
  Similar to other predict methods, which returns predictions from a fitted "\code{bess}" object or a fitted "\code{bess.one}" object.
}
\usage{
  \method{predict}{bess}(object, newdata, type = c("ALL", "AIC", "BIC", "EBIC"),...)
  \method{predict}{bess.one}(object, newdata, ...)
}
\arguments{
  \item{object}{Output from the \code{bess} function or the \code{bess.one} function.}
  \item{newdata}{New data used for prediction.}
  \item{type}{Types of coefficients returned. \code{ type = "AIC"} cooresponds to the predictor with optimal AIC value; \code{ type = "BIC"} cooresponds to the predictor with optimal BIC value;
\code{ type = "EBIC"} cooresponds to the predictor with optimal EBIC value;\code{ type = "ALL"} cooresponds to all predictors in the \code{bess} object. Default is \code{ALL}. }
  \item{\dots}{Additional arguments affecting the predictions produced.}
}
\value{The object returned depends on the types of family.}

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
pred=predict(fit, newdata = data$x)

}
