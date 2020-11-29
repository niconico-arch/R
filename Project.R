if (!require("data.table")) install.packages("data.table")
gc()
rm(list=ls())
library(portfolioBacktest)
library(CVXR)
library(fitHeavyTail)
library(ICSNP)

crypto_file <- "C:/Users/Eileen/SR/R/crypto_data.csv"
index_file <- "C:/Users/Eileen/SR/R/crypto_data_index.csv"
#dat <- read.csv(file)
#dat <- dat%>% select(1,5)
#names(dat)[2]<-"adjusted"

#xts(dat[2], order.by = as.Date(dat$Date,"%d-%B-%y"))

dat_zoo <- read.zoo(crypto_file, index.column = 1, header = T, sep = ",", format = "%m/%d/%Y")
dat_xts <- as.xts(dat_zoo)
index_zoo <- read.zoo(index_file, index.column = 1, header = T, sep = ",", format = "%m/%d/%Y")
index_xts <- as.xts(index_zoo)

dat_xts <- dat_xts[ , colSums(is.na(dat_xts)) == 0]

head(dat_xts)
crypto_data <- list()

crypto_data[['open']] = dat_xts
crypto_data[['high']] = dat_xts
crypto_data[['low']] = dat_xts
crypto_data[['close']] = dat_xts
crypto_data[['volume']] = dat_xts
crypto_data[['adjusted']] = dat_xts
crypto_data[['index']] = index_xts
my_crypto_dataset_list <- stockDataResample(crypto_data, N_sample = 50, T_sample = 365*2, num_datasets = 10)

portfolio_fun <- function(data, lmd = 0.5, alpha = 0.95) {
  X <- as.matrix(diff(log(data$adjusted))[-1])
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
                                     norm1(w)<=1, sum(w) == 1))
  result <- solve(prob)
  return(as.vector(result$getValue(w)))
}



bt <- portfolioBacktest(portfolio_fun, my_crypto_dataset_list, benchmark = c('uniform', 'index'), T_rolling_window = 200, optimize_every = 10, rebalance_every = 1,
                        show_progress_bar = TRUE)
bt_sum <- backtestSummary(bt)
bt_sum$performance_summary
summaryTable(bt_sum, type = "DT", order_col = "Sharpe ratio", order_dir = "desc")
summaryBarPlot(bt_sum, measures = c("Sharpe ratio", "max drawdown"))
summaryTable

data("SP500_symbols")
SP500 <- stockDataDownload(stock_symbols = SP500_symbols, from = "2017-01-01", to = "2020-12-01")
head(SP500)
# resample 10 times from SP500, each with 50 stocks and 2-year consecutive data 
my_dataset_list <- stockDataResample(SP500, N_sample = 50, T_sample = 252*2, num_datasets = 10)

stock_bt <- portfolioBacktest(portfolio_fun, my_dataset_list, benchmark = c('uniform', 'index'), T_rolling_window = 200, optimize_every = 10, rebalance_every = 1,
                        show_progress_bar = TRUE)
stock_bt_sum <- backtestSummary(stock_bt)
stock_bt_sum$performance_summary
summaryTable(stock_bt_sum, type = "DT", order_col = "Sharpe ratio", order_dir = "desc")
summaryBarPlot(stock_bt_sum, measures = c("Sharpe ratio", "max drawdown"))
summaryTable
