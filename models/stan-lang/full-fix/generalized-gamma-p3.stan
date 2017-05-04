data {
  int n_obs;

  vector[n_obs] t_start;
  vector[n_obs] t_stop;
  vector[n_obs] t_truncate;

  real<lower=0> kappa_alpha;
  real<lower=0> kappa_mu;   
  real<lower=0> kappa_beta;
  real<lower=0> sigma_log;    // default is 1/2
  real<lower=0> sigma_logit;    // default is 3/2
}

transformed data {
  vector<lower=0>[n_obs] t_delta;
  t_delta = t_stop-t_start;
}

parameters {
  real alpha_0;
  real mu_0;
  real beta_0;
}


model {
  alpha_0 ~ normal(0,1);
  mu_0 ~ normal(0,1);
  beta_0 ~ normal(0,1);  

  {
  vector[n_obs] alpha;
  vector[n_obs] mu;
  vector[n_obs] beta;
  for ( i in 1:n_obs ) {
    alpha[i] = kappa_alpha * exp(alpha_0/4);
    mu[i] = kappa_mu * exp(mu_0/4);
    beta[i] = kappa_beta * exp(beta_0/4);
  }
 
  for ( i in 1:n_obs ) {
    target += (alpha[i]-.5) * log(alpha[i]) - log(beta[i]) - lgamma(alpha[i]);
    target += sqrt(alpha[i]) * (log(t_delta[i]) - mu[i])/beta[i];
    target += -1 * alpha[i] * 
      exp((log(t_delta[i]) - mu[i])/(beta[i]*sqrt(alpha[i])));
  } 

  }
}




