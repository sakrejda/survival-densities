#include <Rcpp.h>
#include "functions.hpp"

// [[Rcpp::export]]
double gamma_pdf_1S(const double& x, const double& alpha,
  const double& beta 
) {
  double o;
  
  o = functions_model_namespace::gamma_pdf_1S(x, alpha, beta, &Rcpp::Rcout);
  return o;
}

// [[Rcpp::export]]
double gamma_lpdf_1S(const double& x, const double& alpha,
  const double& beta 
) {
  double o;
  
  o = functions_model_namespace::gamma_lpdf_1S(x, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double gamma_cdf_1S(const double& x, const double& alpha,
  const double& beta 
) {
  double o;
  
  o = functions_model_namespace::gamma_cdf_1S(x, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double gamma_lcdf_1S(const double& x, const double& alpha,
  const double& beta 
) {
  double o;
  
  o = functions_model_namespace::gamma_lcdf_1S(x, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double weibull_pdf_1S(const double& x, const double& alpha,
  const double& beta 
) {
  double o;
  
  o = functions_model_namespace::weibull_pdf_1S(x, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double weibull_lpdf_1S(const double& x, const double& alpha,
  const double& beta 
) {
  double o;
  
  o = functions_model_namespace::weibull_lpdf_1S(x, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double weibull_cdf_1S(const double& x, const double& alpha,
  const double& beta 
) {
  double o;
  
  o = functions_model_namespace::weibull_cdf_1S(x, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double weibull_lcdf_1S(const double& x, const double& alpha,
  const double& beta 
) {
  double o;
  
  o = functions_model_namespace::weibull_lcdf_1S(x, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double weibull_uniform_mixture_pdf_1S(const double& x, const double& max,
  const double& q, const double& alpha, const double& beta 
) {
  double o;
  
  o = functions_model_namespace::weibull_uniform_mixture_pdf_1S(x, max, q, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double weibull_uniform_mixture_lpdf_1S(const double& x, const double& max,
  const double& q, const double& alpha, const double& beta 
) {
  double o;
  
  o = functions_model_namespace::weibull_uniform_mixture_lpdf_1S(x, max, q, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double weibull_uniform_mixture_cdf_1S(const double& x, const double& max,
  const double& q, const double& alpha, const double& beta 
) {
  double o;
  
  o = functions_model_namespace::weibull_uniform_mixture_cdf_1S(x, max, q, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double weibull_uniform_mixture_lcdf_1S(const double& x, const double& max, 
  const double& q, const double& alpha, const double& beta 
) {
  double o;
  
  o = functions_model_namespace::weibull_uniform_mixture_lcdf_1S(x, max, q, alpha, beta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double generalized_gamma_pdf_1S(const double& x, const double& alpha,
  const double& beta, const double& nu 
) {
  double o;
  
  o = functions_model_namespace::generalized_gamma_pdf_1S(x, alpha, beta, nu, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double generalized_gamma_lpdf_1S(const double& x, const double& alpha,
  const double& beta, const double& nu 
) {
  double o;
  
  o = functions_model_namespace::generalized_gamma_lpdf_1S(x, alpha, beta, nu, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double generalized_gamma_cdf_1S(const double& x, const double& alpha,
  const double& beta, const double& nu 
) {
  double o;
  
  o = functions_model_namespace::generalized_gamma_cdf_1S(x, alpha, beta, nu, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double generalized_gamma_lcdf_1S(const double& x, const double& alpha,
  const double& beta, const double& nu 
) {
  double o;
  
  o = functions_model_namespace::generalized_gamma_lcdf_1S(x, alpha, beta, nu, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double generalized_gamma_pdf_2S(const double& x, const double& k,
  const double& mu, const double& sigma 
) {
  double o;
  
  o = functions_model_namespace::generalized_gamma_pdf_2S(x, k, mu, sigma, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double generalized_gamma_lpdf_2S(const double& x, const double& k,
  const double& mu, const double& sigma 
) {
  double o;
  
  o = functions_model_namespace::generalized_gamma_lpdf_2S(x, k, mu, sigma, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double generalized_gamma_cdf_2S(const double& x, const double& k,
  const double& mu, const double& sigma 
) {
  double o;
  
  o = functions_model_namespace::generalized_gamma_cdf_2S(x, k, mu, sigma, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double generalized_gamma_lcdf_2S(const double& x, const double& k,
  const double& mu, const double& sigma 
) {
  double o;
  
  o = functions_model_namespace::generalized_gamma_lcdf_2S(x, k, mu, sigma, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double gamma_exp_sum_pdf_1S(const double& x, const double& alpha,
  const double& beta, const double& delta
) {
  double o;
  
  o = functions_model_namespace::gamma_exp_sum_pdf_1S(x, alpha, beta, delta, &Rcpp::Rcout);
  return o;
}

// [[Rcpp::export]]
double gamma_exp_sum_lpdf_1S(const double& x, const double& alpha,
  const double& beta, const double& delta
) {
  double o;
  
  o = functions_model_namespace::gamma_exp_sum_lpdf_1S(x, alpha, beta, delta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double gamma_exp_sum_cdf_1S(const double& x, const double& alpha,
  const double& beta, const double& delta
) {
  double o;
  
  o = functions_model_namespace::gamma_exp_sum_cdf_1S(x, alpha, beta, delta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double gamma_exp_sum_lcdf_1S(const double& x, const double& alpha,
  const double& beta, const double& delta
) {
  double o;
  
  o = functions_model_namespace::gamma_exp_sum_lcdf_1S(x, alpha, beta, delta, &Rcpp::Rcout);
  return o;
}

// [[Rcpp::export]]
double gamma_exp_sum_gamma_mix_pdf_1S(const double& x, const double& q, const double& alpha,
  const double& beta, const double& delta
) {
  double o;
  
  o = functions_model_namespace::gamma_exp_sum_gamma_mix_pdf_1S(x, q, alpha, beta, delta, &Rcpp::Rcout);
  return o;
}

// [[Rcpp::export]]
double gamma_exp_sum_gamma_mix_lpdf_1S(const double& x, const double& q, const double& alpha,
  const double& beta, const double& delta
) {
  double o;
  
  o = functions_model_namespace::gamma_exp_sum_gamma_mix_lpdf_1S(x, q, alpha, beta, delta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double gamma_exp_sum_gamma_mix_cdf_1S(const double& x, const double& q, const double& alpha,
  const double& beta, const double& delta
) {
  double o;
  
  o = functions_model_namespace::gamma_exp_sum_gamma_mix_cdf_1S(x, q, alpha, beta, delta, &Rcpp::Rcout);
  return o;
}


// [[Rcpp::export]]
double gamma_exp_sum_gamma_mix_lcdf_1S(const double& x, const double& q, const double& alpha,
  const double& beta, const double& delta
) {
  double o;
  
  o = functions_model_namespace::gamma_exp_sum_gamma_mix_lcdf_1S(x, q, alpha, beta, delta, &Rcpp::Rcout);
  return o;
}

