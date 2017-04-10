data {
  int n_obs;

  vector[n_obs] t_start;
  vector[n_obs] t_stop;
  vector[n_obs] t_truncate;

  real<lower=0> kappa_alpha;
  real<lower=0> kappa_beta;   
  real<lower=0> kappa_nu;
  real<lower=0> sigma_log;    // default is 1/2
  real<lower=0> sigma_logit;    // default is 3/2
}

transformed data {
  vector<lower=0>[n_obs] t_delta;
  t_delta = t_stop-t_start;
}

parameters {
  real alpha_0;
  real beta_0;
  real nu_0;
}


model {
  alpha_0 ~ normal(0,sigma_log);
  beta_0 ~ normal(0,sigma_log);
  nu_0 ~ normal(0,sigma_logit);  

  {
  vector[n_obs] alpha;
  vector[n_obs] beta;
  vector[n_obs] nu;
  for ( i in 1:n_obs ) {
      alpha[i] = kappa_alpha * exp(alpha_0);
      beta[i] = kappa_beta * exp(beta_0);  
      nu[i] = inv_logit(nu_0);  
  }
  for ( i in 1:n_obs ) {
      target += log(alpha[i]*nu[i]) - log(beta[i]) - lgamma(1/nu[i]) +
        (alpha[i]-1)*(log(t_delta[i]) - log(beta[i])) - 
        pow(t_delta[i]/beta[i],alpha[i]*nu[i]);
      target += log(gamma_p(1/nu[i], pow(t_truncate[i]/beta[i],alpha[i]*nu[i])));
  }
  }
}




