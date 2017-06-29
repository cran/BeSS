\name{bess}
\alias{bess}
\title{Best subset selection}
\description{
  Best subset selection for generalized linear model and Cox's proportional model.
}
\usage{
  bess(x, y, family = c("gaussian", "binomial", "cox"), method = "gsection",
  s.min = 1, s.max, s.list, K.max = 20, max.steps = 15, glm.max = 1e6, cox.max = 20,
  epsilon = 1e-4, normalize = TRUE)
}
\arguments{
  \item{x}{Input matrix,of dimension n x p; each row is an observation vector.}
  \item{y}{Response variable,of length n. For family="binomial" should be a factor with two levels.
           For family="cox", y should be a two-column matrix with columns named 'time' and 'status'.}
  \item{family}{One of the GLM or Cox models. Either "gaussian", "binomial", or "cox", depending on the response.}
  \item{method}{Methods tobe used to select the optimal model size. For method = "\code{sequential}", we solve the best subset selection problem for each \eqn{s} in \eqn{1,2,\dots,s_{max}}. At each model size \eqn{s}, we run the \code{bess} function with a warm start from the last solution with model size \eqn{s-1}. For method = "\code{gsection}", we solve the best subset selection problem with a range non-coninuous model sizes.}
  \item{s.min}{The minimum value of model sizes. Only used for method = "\code{gsection}". Default is 1.}
  \item{s.max}{The maximum value of model sizes. Only used for method = "\code{gsection}". Default is \eqn{\min{p, n/\log(n)}}.}
  \item{s.list}{A list of sequential value representing the model sizes. Only used for method = "\code{sequential}".Default is \eqn{(1,\min{p, n/\log(n)})}.}

  \item{K.max}{The maximum iterations used for method = "\code{gsection}"}
    \item{max.steps}{The maximum number of iterations in \code{bess} function. In linear regression, only a few steps can gurantee the convergence. Default is 15.}
  \item{glm.max}{The maximum number of iterations for solving the maximum likelihood problem on the active set at each step in the primal dual active set algorithm.Only used in the logistic regression for family="binomial". Default is 1e6.}
  \item{cox.max}{The maximum number of iterations for solving the maximum partial likelihood problem on the active set at each step in the primal dual active set algorithm. Only used in Cox's model for family="cox". Default is 20.}

  \item{epsilon}{The tolerance for an early stoping rule in the method "sequential". The early stopping rule is defined as \eqn{\|Y-X\beta\|/n \leq \epsilon}.}
  \item{normalize}{whether to normalize x or not. Default is TRUE.}

}
\details{
  The best subset selection problem with model size \eqn{s} is
  \deqn{\min_\beta -2 logL(\beta) \;\;{\rm s.t.}\;\; \|\beta\|_0 \leq s.}
  In the GLM case, \eqn{logL(\beta)} is the log-likelihood function; In the Cox model, \eqn{logL(\beta)} is the log parital likelihood function.

  For each candiate model size, the best subset selection problem is solved by the primal dual active set(PDAS) algorithm, see Wen et al(2017) for details. This algorithm utilizes an active set updating strategy via primal and dual vairables and fits the sub-model by exploiting the fact that their support set are non-overlap and complementary. For the case of method = "sequential", we run the PDAS algorithm for a list of sequential model sizes and use the estimate from last iteration as a warm start. For the case of method = "\code{gsection}", a golden section search technique is adopted to efficiently determine the optimal model size.
}

\value{
A list with class attribute 'bess' and named components:
  \item{family}{Types of the model: "\code{bess_gaussian}" for linear model,"\code{bess_binomial}" for logistic model and "\code{bess_cox}" for Cox model.}
  \item{beta}{The best fitting coefficients of size \eqn{s=0,1,\dots,p} with the smallest loss function.}
  \item{lambda}{The lambda value in the Lagrangian form of the best subset selection problem with model size of \eqn{s}.}
  \item{deviance}{The value of \eqn{-2\times logL}.}
  \item{nulldeviance}{The value of \eqn{-2\times logL} for null model.}
  \item{AIC}{The value of \eqn{-2\times logL + 2 \|\beta\|_0}.}
  \item{BIC}{The value of \eqn{-2\times logL+ log(n) \|\beta\|_0}.}
  \item{EBIC}{The value of \eqn{-2\times logL+ (log(n)+2\times log(p)) \|\beta\|_0}.}
}


\author{
  Canhong Wen, Aijun Zhang, Shijie Quan, and Xueqin Wang.
}
\references{
 Wen, C., Zhang, A., Quan, S., and Wang, X. (2017) BeSS: A R package for best subset selection in GLM and CoxPH Models, Technical reports.
}

\seealso{
\code{\link{bess.one}}, \code{\link{plot.bess}},
\code{\link{predict.bess}}.
}


\examples{
#--------------linear model--------------#
# Generate simulated data
n <- 500
p <- 20
K <-10
sigma <- 1
rho <- 0.2
data <- gen.data(n, p, family = "gaussian", K, rho, sigma)

# Best subset selection
fit <- bess(data$x, data$y, family = "gaussian")
summary(fit)
coef(fit, sparse=TRUE)  # The estimated coefficients

# Plot solution path and the loss function
plot(fit, type = "both", breaks = TRUE)


#--------------logistic model--------------#

# Generate simulated data
data <- gen.data(n, p, family="binomial", 5, rho, sigma)

# Best subset selection
fit2 <- bess(data$x, data$y, s.list = 1:15, method = "sequential",
             family = "binomial", epsilon = 0)
summary(fit2)
coef(fit2, sparse = TRUE)

# Plot solution path and the loss function
plot(fit2, type = "both", breaks = TRUE, K = 5)

#--------------cox model--------------#

# Generate simulated data
data <- gen.data(n, p, K, rho, sigma, c = 10, family = "cox", scal = 10)

# Best subset selection
fit3 <- bess(data$x, data$y, s.list = 1:15, method = "sequential",
             family = "cox")
coef(fit3, sparse = TRUE)
summary(fit3)

# Plot solution path and the loss function
plot(fit3, type = "both", breaks = TRUE, K = 10)


#----------------------High dimensional linear models--------------------#

p <- 1000
data <- gen.data(n, p, family = "gaussian", K, rho, sigma)

# Best subset selection
fit <- bess(data$x, data$y, method="sequential", family = "gaussian", epsilon = 1e-12)

# Plot solution path
plot(fit, type = "both", breaks = TRUE, K = 10)


}


