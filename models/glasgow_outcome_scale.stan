
functions {
// induced dirichlet prior for cut points from:
// Betancourt, Michael (2019). Ordinal Regression. 
// Retrieved from 
// https://github.com/betanalpha/knitr_case_studies/tree/master/ordinal_regression, 
// commit 1b602ef70a969ff873c04e26101bf44195ce4609.
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}
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
  c ~ induced_dirichlet(rep_vector(1, K), 0);
  // likelihood
  for (n in 1:N) {
    y[n] ~ ordered_logistic(x[n] * beta, c);
  }
}
generated quantities{
  real common_odds_ratio = exp(beta[1]);
  
  // success defined as at least GOS-E 4 "upper severe disability" by investigators
  real control_success = 1 - inv_logit(c[3]);
  real treatment_success = 1 - inv_logit(c[3] - beta[1]);
  real relative_risk = treatment_success / control_success;
  real absolute_risk = treatment_success - control_success;
  
  // posterior predictive distribution
   array[N] int<lower=1, upper=K> y_ppd;
   for (n in 1:N) {
     y_ppd[n] = ordered_logistic_rng(x[n] * beta, c);
     }
}

