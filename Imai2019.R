
## Imai 2019
# install.packages("RcppSMC")

rm(list=ls())

library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(bridgesampling)

data(radiata,package = "RcppSMC")

head(radiata)
plot(radiata)

n <- nrow(radiata)
X <- cbind(rep(1,n),radiata$x1)
y <- matrix(radiata$y,ncol=1)

r_0 <- 0.06
s_0 <- 6
Q_0 <- diag(c(r_0,s_0))
a_0 <- 6
b_0 <- 600^2

## the exact value of marginal likelihood

M <- crossprod(X) + Q_0
R <- diag(n) - X %*% solve(M) %*% t(X)
ML <- pi^(-n/2) * b_0^(a_0/2) * (gamma((n+a_0)/2) / gamma(a_0/2)) *
  (det(Q_0)^(1/2) / det(M)^(1/2)) * as.numeric((t(y)%*%R%*%y+b_0)^(-(n+a_0)/2))
log(ML)
## -312.6449

logML <- -(n/2)*log(pi) + (a_0/2)*log(b_0) + log( gamma((n+a_0)/2)/gamma(a_0/2) ) +
  (1/2)*log( det(Q_0)/det(M) ) - ((n+a_0)/2)*log(t(y)%*%R%*%y+b_0)
logML
## -312.6449

data1 <- list(
  n = n,
  x = c(X[,2]),
  y = c(y),
  mode = 0 ## mode = 1にするとWBIC用 beta = (1/log(n))の事後分布生成
)

model1 <- stan_model("C:/Users/dhojo/Desktop/Imai2019_model1.stan")
fit1 <- sampling(model1, data1, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
stan_rhat(fit1)

data1_wbic <- list(
  n = n,
  x = c(X[,2]),
  y = c(y),
  mode = 1 ## mode = 1にするとWBIC用 beta = (1/log(n))の事後分布生成
)
fit1_wbic <- sampling(model1, data1_wbic, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)


## WBIC
wbic <- function(log_lik){
  wbic <- mean(rowSums(log_lik))
  return(wbic)
}
wbic(rstan::extract(fit1_wbic)$log_lik)
# [1] -308.6954

## WBIC - v_t
## general function variance 汎関数分散
gf_variance <- function(log_lik) sum(colMeans(log_lik^2) - colMeans(log_lik)^2) 

v_t <- (1/(2*log(n)))* gf_variance(rstan::extract(fit1_wbic)$log_lik) ## 特異揺らぎ (singular fluctuation)
wbic(rstan::extract(fit1_wbic)$log_lik)-v_t ## WBIC-v_t

### bridge sampling

logML_bs1 <- bridge_sampler(fit1,method = "normal",repetitions = 10, cores = 4, silent = TRUE)
logML_bs1
# Median of 10 bridge sampling estimates
# of the log marginal likelihood: -310.5496
# Range of estimates: -310.5525 to -310.5442
# Interquartile range: 0.00474
# Method: normal

logML_warp1 <- bridge_sampler(fit1,method = "warp3",repetitions = 10, cores = 4, silent = TRUE)
logML_warp1

# Median of 10 bridge sampling estimates
# of the log marginal likelihood: -310.5488
# Range of estimates: -310.5501 to -310.5458
# Interquartile range: 0.00158
# Method: warp3


### Thermodynamic integration
model_ti <- stan_model("C:/Users/dhojo/Desktop/Imai2019_model1_temp.stan")

K = 11
temperature <- rep(NA,K)
for(j in 1:K) {
  temperature[j] <- ((j-1)/(K-1))^(1/.3)
}
plot(temperature)
temperature

data1_0.0 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[1]) ## prior only
data1_0.1 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[2])
data1_0.2 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[3])
data1_0.3 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[4])
data1_0.4 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[5])
data1_0.5 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[6])
data1_0.6 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[7])
data1_0.7 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[8])
data1_0.8 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[9])
data1_0.9 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[10])
data1_1.0 <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[11]) ## posterior

fit_0.0 <- sampling(model_ti, data1_0.0, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
fit_0.1 <- sampling(model_ti, data1_0.1, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
fit_0.2 <- sampling(model_ti, data1_0.2, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
fit_0.3 <- sampling(model_ti, data1_0.3, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
fit_0.4 <- sampling(model_ti, data1_0.4, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
fit_0.5 <- sampling(model_ti, data1_0.5, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
fit_0.6 <- sampling(model_ti, data1_0.6, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
fit_0.7 <- sampling(model_ti, data1_0.7, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
fit_0.8 <- sampling(model_ti, data1_0.8, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
fit_0.9 <- sampling(model_ti, data1_0.9, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
fit_1.0 <- sampling(model_ti, data1_1.0, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)

## mean(rowSums())するのでwbic関数使っています｡
ll_0.0 <- wbic(rstan::extract(fit_0.0)$log_lik)
ll_0.1 <- wbic(rstan::extract(fit_0.1)$log_lik)
ll_0.2 <- wbic(rstan::extract(fit_0.2)$log_lik)
ll_0.3 <- wbic(rstan::extract(fit_0.3)$log_lik)
ll_0.4 <- wbic(rstan::extract(fit_0.4)$log_lik)
ll_0.5 <- wbic(rstan::extract(fit_0.5)$log_lik)
ll_0.6 <- wbic(rstan::extract(fit_0.6)$log_lik)
ll_0.7 <- wbic(rstan::extract(fit_0.7)$log_lik)
ll_0.8 <- wbic(rstan::extract(fit_0.8)$log_lik)
ll_0.9 <- wbic(rstan::extract(fit_0.9)$log_lik)
ll_1.0 <- wbic(rstan::extract(fit_1.0)$log_lik)

ll <- data.frame(logml = c(ll_0.0,ll_0.1,ll_0.2,ll_0.3,ll_0.4,ll_0.5,ll_0.6,ll_0.7,ll_0.8,ll_0.9,ll_1.0),
                 temperature = temperature)
ll %>% 
  ggplot(aes(x=temperature,y=ml)) + 
  geom_point() +
  theme_bw()


logML_ti <- c(rep(NA,K-1))
for(k in 2:K){
  logML_ti[k-1] = ((temperature[k] - temperature[k-1])/2) * (ll[ll$temperature == temperature[k],"logml"]+ll[ll$temperature == temperature[k-1],"logml"])
}

sum(logML_ti) ## log ML by TI method
# [1] -311.1424

## sim

curve(-log(x),col="red")
# for (i in seq(0.1,1.0,0.1)) {
#   curve(dbeta(x,i,1))
#   par(new=TRUE)
# }





