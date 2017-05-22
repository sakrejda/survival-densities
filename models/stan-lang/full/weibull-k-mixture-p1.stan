data {
  int n_obs;
  int n_levels;
  int n_groups[n_levels];

  vector[n_obs] t_start;
  vector[n_obs] t_stop;
  vector[n_obs] t_truncate;

  int n_components;
  real<lower=0> kappa_alpha;
  real<lower=0> kappa_beta;   
  real<lower=0> sigma_log;    // default is 1/2
  real<lower=0> sigma_logit;  // default is 3/2
}

transformed data {
  vector[n_obs] t_delta;
  t_delta = t_stop-t_start;
}

parameters {
  vector[n_components] alpha_0;
  ordered[n_components] beta_0;
  vector[n_components] q_0;

}

transformed parameters {
}

model {
  alpha_0 ~ normal(0,sigma_log);
  beta_0 ~ normal(0,sigma_log);
  q_0 ~ normal(0,sigma_logit);  


  {
  matrix[n_components, n_obs] alpha;
  matrix[n_components, n_obs] beta;
  matrix[n_components, n_obs] q_lin;
  matrix[n_components, n_obs] q;
  vector[n_components] mix;
  for ( i in 1:n_obs ) {
    for ( k in 1:n_components ) {
      alpha[k,i] = kappa_alpha * exp(alpha_0[k]);
      beta[k,i] = kappa_beta * exp(beta_0[k]);  
      q_lin[k,i] = q_0[k];
    }
  }
  for ( i in 1:n_obs ) {
    q[,i] = softmax(q_lin[,i]);
  }
  for ( i in 1:n_obs ) {
    for ( k in 1:n_components ) {
      mix[k] = log(q[k,i]) + 
        weibull_lpdf(t_delta[i] | alpha[k,i], beta[k,i]) - 
        weibull_lcdf(t_truncate[i] | alpha[k,i], beta[k,i]);
    }
    target += log_sum_exp(mix);
  }
  }
}




