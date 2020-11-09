install.packages('portfolioBacktest')

rm(list=ls())
gc()
library(portfolioBacktest)
library(PerformanceAnalytics)
library(CVXR)

data("dataset10")

MVP <- function(mu, Sigma, upper, lmd) {
  w <- Variable(nrow(Sigma))
  prob <- Problem(Maximize(t(mu) %*% w - lmd*quad_form(w, Sigma)),
                  constraints = list(w >= 0, w <= upper, sum(w) == 1))
  result <- CVXR::solve(prob)
  w <- as.vector(result$getValue(w))
  names(w) <- colnames(Sigma)
  return(w)
}


l = rep(0, 50)
S = rep(NA, 50)
result <- list()

for (d in dataset10[1:5]){
  X_log <- diff(log(d$adjusted))[-1]
  N <- ncol(X_log)
  mu <- colMeans(X_log)
  Sigma <- cov(X_log)
  u = rep(1, 50)
  for (i in 1:length(u)){
    w_ <- MVP(mu, Sigma, u[i], l[i])
    ret <- xts(X_log %*% as.numeric(w_), index(X_log))
    S[i] <- SharpeRatio(R=ret, FUN = "StdDev")
  }
  result <- c(result, cbind(l, u, S))
}
