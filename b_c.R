#---------------------------------- CASE B ---------------------------------------------(1)
# input: symbols: 20 symbols from Tickers.csv equally weighted portfolio
# output: mean return and volatility using covariance matrix

library(quantmod)
library(xts)
library(PerformanceAnalytics)
library(PortfolioAnalytics)
library(tseries)

library(DEoptim)
library(ROI)
require("zoo")
#require(ROI.plugin.glpk)
#require(ROI.plugin.quadprog)

# all 20 symbols from Ticker list
symbols <- c("MU", "CMCSA", "MSFT", "QQQ", "AAPL", "NVDA", "INTC", "PYPL", 
             "CSCO", "TVIX", "FB", "SIRI", "AMAT", "QCOM", "XNET", "FOXA", 
             "JD", "AABA", "DISCA", "SBUX")

# length is 20
l <- length(symbols)

# create vector of 5% weights for all symbols
weights <- rep(1/l, l)

# all 20 symbols prices in one place, first initializing prices variable
prices <- xts()

# loop through all 20 symbols and get their prices for 2016
for (i in 1:l) {
  symbol <- symbols[i]
  data <- getSymbols(symbol)
  prices <- merge(prices, Ad(get(symbols[i]))) # Op, Cl, etc.
  rm(symbol)
}

# calculate returns for each asset using Return.calculate()  
returns <- Return.calculate(prices['2016'])
# Set names
colnames(returns) <- symbols

# Remove the first row of returns - [NA values]
returns <- returns[-1, ]

#---------------------------------------------------------------------------------------(2)
# Covariance matrix
sigma <- cov(returns)

# Create a weights matrix w
w <- as.matrix(weights)

# sum of returns for each ticker in 2016
total_returns <- as.matrix(colSums(returns))

# total portfolio's return in 2016:
# checking by sum(Return.portfolio(R = returns, weights = weights, geometric = FALSE))
total <- t(w) %*% total_returns
total

# variance of portfolio return is:
port_var <- t(w) %*% sigma %*% w
port_var

#---------------------------------- CASE C ---------------------------------------------(3)

# input: symbols: 20 symbols from Tickers.csv equally weighted portfolio
# output: optimal w1...w20 weights that minize the risk by using optimization method ROI,  
# build efficient frontier
# goal: Minimize variance with ROI

# 1. Create portfolio specification
port_spec <- portfolio.spec(assets = colnames(returns))

# 2. Add constraints
# Add a full investment constraint, sum of weights equal to 1
port_spec <- add.constraint(portfolio = port_spec, type = "full_investment")

# Add a long only constraint, weight of an asset is between 0 and 1
port_spec <- add.constraint(portfolio = port_spec, type = "long_only")

# 3. Add an objective to minimize portfolio variance
port_spec <- add.objective(portfolio = port_spec, type = "risk", name = "var")

# 4. Solve the optimization problem with ROI
opt <- optimize.portfolio(R = returns, portfolio = port_spec, optimize_method = "ROI")

# Extract the optimal weights
opt_w <- c(extractWeights(opt))
opt_w
sum(opt_w)

total_optim <- t(as.matrix(opt_w)) %*% total_returns

#---------------------------------------------------------------------------------------(4)

# Analyzing performances of two portfolio cases: equally weighted and weight-optimized

case_b <- Return.portfolio(R = returns, weights = weights, geometric = FALSE)
case_c <- Return.portfolio(R = returns, weights = opt_w, geometric = FALSE)

# merge case b and c in to one
ps = merge(case_b, case_c)
names(ps) = c("case B", "case C")

# chart comaprison of two cases
charts.PerformanceSummary(ps, colorset=c("red", "blue"), legend.loc="topleft")

#---------------------------------------------------------------------------------------(5)

# Computing the efficient frontier using a grid of target returns

# Calculate each stocks mean returns
stockmu <- colMeans(returns)

# Create a grid of target values
grid <- seq(from = 0.01, to = max(stockmu), length.out = 50)

# Create empty vectors to store means and deviations
vpm <- vpsd <- rep(NA, length(grid))

# Create an empty matrix to store weights
mweights <- matrix(NA, 50, 20)

# Create your for loop
for(i in 1:length(grid)) {
  optimum <- portfolio.optim(x = returns, pm = grid[i])
  vpm[i] <- optimum$pm
  vpsd[i] <- optimum$ps
  mweights[i, ] <- optimum$pw
}

plot(vpsd, vpm)



