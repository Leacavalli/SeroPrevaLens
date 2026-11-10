// Mixture1D_priors.stan
// Mixture1D.stan with the prior hyperparameters and the parameter limits passed as data,
// matching the prior-sensitivity page. With the Table 1 values of O'Driscoll et al. (2025)
// it reproduces Mixture1D.stan exactly. Limits may be infinite (pass Inf / -Inf).
data {
  int<lower=1> N;              // number of observations
  vector[N] y;                 // log(RFI) or log(MFI)

  real<lower=0> pi_a;  real<lower=0> pi_b;   // pi ~ beta(pi_a, pi_b)
  real mu0_loc;   real<lower=0> mu0_scale;   // normal(loc, scale), truncated to the limits
  real muR_loc;   real<lower=0> muR_scale;
  real sd0_loc;   real<lower=0> sd0_scale;
  real sd1_loc;   real<lower=0> sd1_scale;

  // parameter limits
  real mu0_lower;  real<lower=mu0_lower> mu0_upper;
  real<lower=0> muR_lower;  real<lower=muR_lower> muR_upper;
  real<lower=0> sd0_lower;  real<lower=sd0_lower> sd0_upper;
  real<lower=0> sd1_lower;  real<lower=sd1_lower> sd1_upper;
}
parameters {
  real<lower=0, upper=1>                 pi;   // prevalence
  real<lower=mu0_lower, upper=mu0_upper> mu0;  // mean of negative titers
  real<lower=muR_lower, upper=muR_upper> muR;  // mean increase in positive titers
  real<lower=sd0_lower, upper=sd0_upper> sd0;  // SD of negative titers
  real<lower=sd1_lower, upper=sd1_upper> sd1;  // SD of positive titers
}
transformed parameters {
  real mu1 = mu0 + muR;        // mean of positive titers
}
model {
  pi  ~ beta(pi_a, pi_b);
  mu0 ~ normal(mu0_loc, mu0_scale);
  muR ~ normal(muR_loc, muR_scale);
  sd0 ~ normal(sd0_loc, sd0_scale);
  sd1 ~ normal(sd1_loc, sd1_scale);

  // Likelihood: (1 - pi) N(mu0, sd0) + pi N(mu1, sd1)
  for (n in 1:N) {
    target += log_mix(pi,
                      normal_lpdf(y[n] | mu1, sd1),   // positive, weight pi
                      normal_lpdf(y[n] | mu0, sd0));  // negative, weight 1 - pi
  }
}
generated quantities {
  vector[N] log_lik;
  vector[N] p_pos;
  real sumloglik;
  real cutoff = (mu0 + mu1) / 2;

  for (n in 1:N) {
    real lp1 = log(pi)   + normal_lpdf(y[n] | mu1, sd1);
    real lp0 = log1m(pi) + normal_lpdf(y[n] | mu0, sd0);
    log_lik[n] = log_sum_exp(lp1, lp0);
    p_pos[n]   = exp(lp1 - log_lik[n]);
  }
  sumloglik = sum(log_lik);
}
