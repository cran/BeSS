## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = F,
  message = F
)

## ---- warning=F, message=FALSE, eval=TRUE-------------------------------------
library(BeSS)
data("trim32", package = "BeSS")

## -----------------------------------------------------------------------------
y <- trim32[, 1]
x <- as.matrix(trim32[, -1])
lm.bsrr <- bess(x, y, type = "bsrr")

## ---- eval=F------------------------------------------------------------------
#  coef(lm.bsrr, sparse = TRUE)

## ---- eval=F------------------------------------------------------------------
#  predict.bsrr <- predict(lm.bsrr, newx = x)

## ---- warning=FALSE, message = FALSE------------------------------------------
data("duke")
y <- duke$y
x <- as.matrix(duke[, -1])

## -----------------------------------------------------------------------------
logi.bsrr <- bess(x, y, family = "binomial", type = "bsrr")

## -----------------------------------------------------------------------------
plot(logi.bsrr)

## ---- warning=FALSE, message = FALSE------------------------------------------
data(LymphomaData, package = "HCmodelSets")

x <- t(patient.data$x)
y <- patient.data$time
status <- patient.data$status

## -----------------------------------------------------------------------------
cox.bsrr <- bess(x, cbind(y, status), family = "cox", type = "bsrr")

## -----------------------------------------------------------------------------
summary(cox.bsrr)

## ---- eval=F------------------------------------------------------------------
#  lm.bsrr.ebic <- bess(x, y, type = "bsrr", tune = "ebic")

## ---- eval=F------------------------------------------------------------------
#  lm.bsrr.cv <- bess(x, y, type = "bsrr", tune = "cv", nfolds = 5)

## ---- eval=F------------------------------------------------------------------
#  my.lambda.list <- exp(seq(log(10), log(0.01), length.out = 10))
#  my.s.list <- 1:10
#  
#  lm.bsrr.seq <- bess(x, y, type = "bsrr", method = "sequential", s.list = my.s.list,
#                     lambda.list = my.lambda.list)

## ---- eval=F------------------------------------------------------------------
#  my.s.min <- 1
#  my.s.max <- 10
#  my.lambda.min <- 0.01
#  my.lambda.max <- 10
#  
#  lm.bsrr.powell <- bess(x, y, type = "bsrr", method = "pgsection",
#                        s.min = my.s.min, s.max = my.s.max,
#                        lambda.min = my.lambda.min, lambda.max = my.lambda.max)

## ---- eval=F------------------------------------------------------------------
#  lm.bsrr.screening <- bess(x, y, type = "bsrr", screening.num = round(nrow(x)/log(nrow(x))))

