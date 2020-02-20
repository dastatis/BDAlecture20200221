data{
  int<lower = 0> n;
  vector[n] x;
  vector[n] y;
  int<lower=0, upper=1> mode; // 0 = normal sampling, 1 = WBIC sampling
}

transformed data {
  // https://tbz533.blogspot.com/2017/04/how-to-calculate-wbic-with-stan.html
    real<lower=0, upper=1> watanabe_beta; // WBIC parameter
    if ( mode == 0 ) {
        watanabe_beta = 1.0;
    }
    else { // mode == 1
        watanabe_beta = 1.0/log(n);
    }
}

parameters {
  real alpha;
  real beta;
  real<lower = 0> tau;
}

transformed parameters {
  vector[n] mu;
  for (i in 1:n){
    mu[i] = alpha + beta * (x[i] - mean(x));
  }
}

model {
  target += gamma_lpdf(tau | 6, 600^2);
  target += normal_lpdf(alpha | 3000, sqrt(tau * 0.06)^-1);
  target += normal_lpdf(beta | 185, sqrt(tau * 6)^-1);
  target += watanabe_beta * normal_lpdf(y|mu, sqrt(tau)^-1);
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n){
    log_lik[i] = normal_lpdf(y[i]|mu[i], sqrt(tau)^-1);
  }
}
