// One-factor confirmatory factor analysis (marginal likelihood)
// Integrates out the latent factor analytically
data {
  int<lower=1> N;         // number of observations
  int<lower=1> P;         // number of indicators
  matrix[N, P] Y;         // response matrix (N x P)
}

parameters {
  vector[P] alpha;        // intercepts
  vector[P] lambda;       // factor loadings (unconstrained, allows sign switching)
  vector<lower=0>[P] sigma;  // residual standard deviations
}

transformed parameters {
  // Marginal covariance: V = lambda * lambda' + diag(sigma^2)
  matrix[P, P] V;
  V = tcrossprod(to_matrix(lambda, P, 1)) + diag_matrix(square(sigma));
}

model {
  // Priors
  alpha ~ normal(0, 5);
  lambda ~ normal(0, 1);     // symmetric prior: allows sign switching
  sigma ~ normal(0, 2);

  // Marginal likelihood (integrated over latent factor)
  for (j in 1:N) {
    Y[j] ~ multi_normal(alpha, V);
  }
}

generated quantities {
  // Pointwise MARGINAL log-likelihoods for DIC_i, WAIC, LOO
  vector[N] log_lik;
  for (j in 1:N) {
    log_lik[j] = multi_normal_lpdf(Y[j] | alpha, V);
  }
}
