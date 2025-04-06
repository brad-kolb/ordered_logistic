data {
  int<lower=2> K; // number of ordinal categories
  int<lower=0> N; //  number of observations
  int<lower=1> D; // number of coefficients
  array[N] int<lower=1, upper=K> y; // observed ordinal outcomes
  array[N] row_vector[D] x; // predictors
}
parameters {
  vector[D] beta; // latent effect
  ordered[K - 1] c; // internal cut points
}
model {
  // prior
  beta ~ normal(0, 1);
  // likelihood
  for (n in 1:N) {
    y[n] ~ ordered_logistic(x[n] * beta, c);
  }
}
generated quantities{
  real common_odds_ratio = exp(beta[1]);
  
  real control_mrs02 = 1 - inv_logit(c[4]);
  real treatment_mrs02 = 1 - inv_logit(c[4] - beta[1]);
  real rr_mrs02 = treatment_mrs02 / control_mrs02;
  
  real control_mrs03 = 1 - inv_logit(c[3]);
  real treatment_mrs03 = 1 - inv_logit(c[3] - beta[1]);
  real rr_mrs03 = treatment_mrs03 / control_mrs03;
  
  // prior predictive check
  // array[N] int<lower=1, upper=K> y_rep;
  // for (n in 1:N) {
  //  y_rep[n] = ordered_logistic_rng(x[n] * beta, c);
  // }
}


