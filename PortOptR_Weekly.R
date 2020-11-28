gc()
rm(list=ls())
library(portfolioBacktest)
library(CVXR)
library(fitHeavyTail)
library(ICSNP)
#if (!require("ICSNP")) install.packages("ICSNP")
#xfun::session_info('DT')

data("dataset10")

data("SP500_symbols")
SP500 <- stockDataDownload(stock_symbols = SP500_symbols, from = "2017-12-01", to = "2018-12-01")
# resample 10 times from SP500, each with 50 stocks and 2-year consecutive data 
my_dataset_list <- stockDataResample(SP500, N_sample = 50, T_sample = 252*2, num_datasets = 10)


# names(SP500)
# names(my_dataset_list$"dataset 1")
# my_dataset_list[1]
 # GMVP <- function(Sigma) {
 #   w <- Variable(nrow(Sigma))
 #   prob <- Problem(Minimize(quad_form(w, Sigma)), 
 #                   constraints = list(norm1(w)<=1, sum(w) == 1))
 #   result <- CVXR::solve(prob)
 #   w <- as.vector(result$getValue(w))
 #   names(w) <- colnames(Sigma)
 #   return(w)
 # }
 # 
 # MDCP <- function(Sigma) {
#    C <- diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma)))
#    colnames(C) <- colnames(Sigma)
#    return(GMVP(Sigma = C))
#  }
# 
# # Tyler estimator with shrinkage
# estimateTylerShrinkage <- function(X, R_target, alpha, verbose = FALSE) {
#   max_iter <- 100
#   error_th_Sigma <- 1e-3
#   
#   #Gaussian initial point
#   Sigma <- 1/(1+alpha)*cov(X) + alpha/(1+alpha)*R_target
#   Sigma <- Sigma/sum(diag(Sigma))
#   #loop
#   obj_value_record <- Sigma_diff_record <- rep(NA, max_iter)
#   for (k in 1:max_iter) {
#     Sigma_prev <- Sigma
#     
#     #Tyler update
#     weights <- 1/rowSums(X * (X %*% solve(Sigma)))
#     obj_value_record[k] <- - (N/2)*sum(log(weights)) + (T/2)*sum(log(eigen(Sigma)$values))
#     Sigma <- (N/T) * crossprod( sqrt(weights)*X )
#     Sigma <- 1/(1+alpha)*Sigma + alpha/(1+alpha)*R_target
#     Sigma <- Sigma/sum(diag(Sigma))
#     
#     #stopping criterion
#     Sigma_diff_record[k] <- norm(Sigma - Sigma_prev, "F")/norm(Sigma_prev, "F")
#     if (Sigma_diff_record[k] < error_th_Sigma)
#       break
#   }
#   obj_value_record <- obj_value_record[1:k]
#   Sigma_diff_record <- Sigma_diff_record[1:k]
#   if (verbose)
#     plot(obj_value_record, type = "l", col = "blue",
#          main = "Convergence of objective value", xlab = "iterations", ylab = "obj value")
#   
#   #recover missing scaling factor
#   sigma2 <- apply(X, 2, var)
#   d <- diag(Sigma)
#   kappa <- sum(sigma2*d)/sum(d*d)
#   Sigma <- kappa*Sigma
#   
#   return(Sigma)
# }
# 
# portfolio_fun <- function(data) {
#   X <- as.matrix(diff(log(data$adjusted))[-1])  # compute log returns
#   N <- ncol(X) # number of stocks
#   T <- nrow(X) #number of days
#   
#   #estimation of R using Tyler estimator with shrinkage
#   mu_gmedian <- ICSNP::spatial.median(X)
#   X_ <- X - rep(mu_gmedian, each = T)
#   alpha <- min(0.01*N/T, 1)
#   R_target <- mean(diag(cov(X))) * diag(N)
#   R_Tyler_shrinked <- estimateTylerShrinkage(X_, R_target, alpha)
# 
#   #Optimize with Maximum decorrelation portfolio with shrinkage estimated Sigma
#   w <- MDCP(R_Tyler_shrinked)
#   return(w)
# }

#mean CVaR portfolio
portfolio_fun <- function(data, lmd = 0.5, alpha = 0.95) {
  X <- as.matrix(diff(log(data$adjusted))[-1])
  #X <- as.matrix(diff(log(dataset10$`dataset 1`$adjusted))[-1])
  N <- ncol(X) # number of stocks
  T <- nrow(X) #number of days
  mu <- colMeans(X)
  
  # variables
  w <- Variable(N)
  z <- Variable(T)
  zeta <- Variable(1)
  
  # problem
  prob <- Problem(Maximize(t(w) %*% mu - lmd*zeta - (lmd/(T*(1-alpha))) * sum(z)),
                  constraints = list(z >= 0, z >= -X %*% w - zeta,
                                     norm1(w)<=1, sum(w) == 1, w>=0))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
  # prob <- Problem(Maximize(t(w) %*% mu - lmd*zeta - (lmd/(T*(1-alpha)))*sum(z)), constraints = list(z>=0, z>=-X%*%w-zeta, norm1(w)<=1, sum(w) == 1))

}

bt <- portfolioBacktest(portfolio_fun, my_dataset_list, benchmark = c('uniform', 'index'), T_rolling_window = 200, optimize_every = 10, rebalance_every = 1,
                        show_progress_bar = TRUE)
bt_sum <- backtestSummary(bt)
bt_sum$performance_summary
summaryTable(bt_sum, type = "DT", order_col = "Sharpe ratio", order_dir = "desc")
summaryBarPlot(bt_sum, measures = c("Sharpe ratio", "max drawdown"))
summaryTable

