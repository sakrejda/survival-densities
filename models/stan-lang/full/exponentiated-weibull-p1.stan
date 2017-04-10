data {
  int n_obs;

  vector[n_obs] t_start;
  vector[n_obs] t_stop;
  vector[n_obs] t_truncate;

  real<lower=0> kappa_alpha;
  real<lower=0> kappa_theta;   
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
  real theta_0;
  real beta_0;
}


model {
  alpha_0 ~ normal(1,sigma_log);
  theta_0 ~ normal(0,sigma_log);
  beta_0 ~ normal(1,sigma_log);  

  {
  vector[n_obs] alpha;
  vector[n_obs] theta;
  vector[n_obs] beta;
  real tba;
  for ( i in 1:n_obs ) {
      alpha[i] = kappa_alpha * exp(alpha_0);
      theta[i] = kappa_theta * exp(theta_0);  
      beta[i] = kappa_beta * exp(beta_0);  
  }
  for ( i in 1:n_obs ) {
    tba = (t_delta[i]/beta[i])^alpha[i];
    target += log(alpha[i]) + log(theta[i]) - log(beta[i]);
    target += -tba;
    target += (alpha[i] - 1) * (log(t_delta[i])-log(beta[i]));
    target += (theta[i] - 1) * (1-exp(-tba));
  }
  }
}





