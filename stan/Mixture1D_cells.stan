// Mixture1D_cells.stan
// Mixture1D_priors.stan with one prevalence per cell (e.g. subdistrict x age group);
// mu0, muR, sd0 and sd1 are shared by all cells. With G = 1 it is Mixture1D_priors.stan.
data {
  int<lower=1> N;                          // number of observations
  vector[N] y;                             // log(RFI) or log(MFI)
  int<lower=1> G;                          // number of cells
  array[N] int<lower=1, upper=G> cell;     // cell of each observation

  real<lower=0> pi_a;  real<lower=0> pi_b;   // pi[g] ~ beta(pi_a, pi_b)
  real mu0_loc;   real<lower=0> mu0_scale;
  real muR_loc;   real<lower=0> muR_scale;
  real sd0_loc;   real<lower=0> sd0_scale;
  real sd1_loc;   real<lower=0> sd1_scale;

  real mu0_lower;  real<lower=mu0_lower> mu0_upper;
  real<lower=0> muR_lower;  real<lower=muR_lower> muR_upper;
  real<lower=0> sd0_lower;  real<lower=sd0_lower> sd0_upper;
  real<lower=0> sd1_lower;  real<lower=sd1_lower> sd1_upper;
}
parameters {
  vector<lower=0, upper=1>[G]            pi;   // prevalence in each cell
  real<lower=mu0_lower, upper=mu0_upper> mu0;
  real<lower=muR_lower, upper=muR_upper> muR;
  real<lower=sd0_lower, upper=sd0_upper> sd0;
  real<lower=sd1_lower, upper=sd1_upper> sd1;
}
transformed parameters {
  real mu1 = mu0 + muR;
}
model {
  pi  ~ beta(pi_a, pi_b);
  mu0 ~ normal(mu0_loc, mu0_scale);
  muR ~ normal(muR_loc, muR_scale);
  sd0 ~ normal(sd0_loc, sd0_scale);
  sd1 ~ normal(sd1_loc, sd1_scale);
  for (n in 1:N)
    target += log_mix(pi[cell[n]], normal_lpdf(y[n] | mu1, sd1), normal_lpdf(y[n] | mu0, sd0));
}
generated quantities {
  vector[N] log_lik;                       // for loo::loo()
  for (n in 1:N)
    log_lik[n] = log_mix(pi[cell[n]], normal_lpdf(y[n] | mu1, sd1), normal_lpdf(y[n] | mu0, sd0));
}
