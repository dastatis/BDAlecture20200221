
## Imai 2019
# install.packages("RcppSMC")


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

model1 <- stan_model("C:/Users/dhojo/Desktop/imai/Imai2019_model1.stan")
fit1 <- sampling(model1, data1, iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
stan_rhat(fit1)

waic_fit1 <- loo::waic(rstan::extract(fit1)$log_lik)
waic_fit1$pointwise ## elpd(expected log predictive density) = lpd() - p_waic(汎関数分散)

apply(rstan::extract(fit1)$log_lik,2,mean) ## データポイント毎 lpd
apply(rstan::extract(fit1)$log_lik,2,var) ## データポイント毎 p_waic

# ## general function variance 汎関数分散 = p_waic
gf_variance <- function(log_lik) sum(colMeans(log_lik^2) - colMeans(log_lik)^2)
gf_variance(rstan::extract(fit1)$log_lik) ## p_waic
sum(apply(rstan::extract(fit1)$log_lik,2,mean)) ## lpd


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
## WBIC
wbic <- function(log_lik){
  wbic <- mean(rowSums(log_lik))
  return(wbic)
}

model_ti <- stan_model("C:/Users/dhojo/Desktop/imai/Imai2019_model1_temp.stan")

K = 6
temperature <-list_along(K)
for(j in 1:K) {
  temperature[[j]] <- ((j-1)/(K-1))^(1/.3)
}
plot(unlist(temperature))
unlist(temperature)

data_ti <- list_along(K)
fit_ti <- list_along(K)
logml_ti_temp <- list_along(K)
for(j in 1:K){
  print(j)
  data_ti[[j]] <- list(n = n, x = c(X[,2]), y = c(y), temperature = temperature[[j]])
  fit_ti[[j]] <- sampling(model_ti, data_ti[[j]], iter = 2000, warmup = 1000, thin  = 1, chains = 4 , seed = 1234)
  logml_ti_temp[[j]] <- wbic(rstan::extract(fit_ti[[j]])$log_lik)
}

ll <- data.frame(logml = unlist(logml_ti_temp),
                 temperature = unlist(temperature))
ll %>% 
  ggplot(aes(x=temperature,y=logml)) + 
  geom_point() +
  theme_bw()


logML_ti <- c(rep(NA,K-1))
for(k in 2:K){
  logML_ti[k-1] = ((temperature[[k]] - temperature[[k-1]])/2) * (ll[ll$temperature == temperature[[k]],"logml"]+ll[ll$temperature == temperature[[k-1]],"logml"])
}

sum(logML_ti) ## log ML by TI method
# [1] -312.7311 K = 5
# [1] -311.1424 K = 11
# [1] -310.795  K = 16
# [1] -310.6975 K = 21
# [1] -310.6057 K = 51
# [1] -310.5809 K = 101



