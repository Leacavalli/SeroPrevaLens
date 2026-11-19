// multiSero_noStrata_v2_antigenLoc_seroVec.stan
// =================================================================
// Multivariate serological mixture model – no strata
//
// v2_antigenLoc changes over multiSero_noStrata_improved_v2.stan:
//   - Per-ANTIGEN prior LOCATIONS for mu0/sd0 (vector[nP]) and
//     per-PRESENT-PATHOGEN prior locations for mu1/sd1 (vector[nPp]),
//     instead of a single pooled scalar location shared by all antigens.
//
//     Motivation: the pooled-median prior forces every antigen's
//     negative component toward one common baseline (mu0 ~ N(median, .),
//     sd0 ~ N(median, .)). That is fine when all antigens share a
//     baseline, but wrong for outlier antigens whose negative population
//     sits elsewhere (e.g. vaccine-shifted YFV NS1: mu_neg ~5.1, sd_neg
//     ~1.2, far above the pooled median). Centring each antigen's prior
//     on its OWN empirical GMM value removes that mis-specification and
//     restores identifiability of the outlier's sero/sd1 without loosening
//     the well-behaved antigens.
//
//   - ONLY the prior-location TYPES change (scalar -> vector). The model
//     block is byte-identical to _v2: `vec ~ normal(vec_loc, real_scale)`
//     vectorises element-wise. Scales remain scalar (weakly informative,
//     shared).
//
// =================================================================

data {
  int<lower=1> N;                   // individuals
  int<lower=1> nP;                  // total antigens
  int<lower=1> nPp;                 // present pathogens
  array[nP] int pres;              // 1 if present pathogen, else 0
  array[N] vector[nP] y;           // log-transformed antibody titers
  int<lower=1> nC;                 // infection-status combinations (2^nPp)
  matrix[nC, nP] infM;            // infection-status indicator per combo
  array[nC] int npos;             // # positive pathogens per combo
  array[nC, nP] int wpos;         // indices of positive antigens (0-padded)
  array[nC, nP] int wneg;         // indices of negative antigens (0-padded)

  // Configurable prior hyperparameters
  // ---- LOCATIONS are now per-antigen / per-present-pathogen VECTORS ----
  vector[nP]  mu0_prior_loc;              // was: real
  real<lower=0> mu0_prior_scale;
  vector[nPp] mu1_prior_loc;              // was: real
  real<lower=0> mu1_prior_scale;
  vector<lower=0>[nP]  sd0_prior_loc;     // was: real<lower=0>
  real<lower=0> sd0_prior_scale;
  vector<lower=0>[nPp] sd1_prior_loc;     // was: real<lower=0>
  real<lower=0> sd1_prior_scale;
  // ---------------------------------------------------------------------
  vector<lower=0>[nPp] sero_alpha;   // was: real  (per-present Beta prior; fill rep(x,nPp) to recover the shared scalar prior)
  vector<lower=0>[nPp] sero_beta;    // was: real
  real<lower=0> phi_rate;
  real<lower=0> rho_alpha;
  real<lower=0> rho_beta;
}

transformed data {
  // Pre-compute data matrix for vectorised likelihood
  matrix[nP, N] Y;
  for (n in 1:N) Y[ : , n] = y[n];
  real log_2pi = log(2 * pi());

  // Pre-compute complements of infM for log_theta calculation
  matrix[nC, nP] infM_inv = 1.0 - infM;

  // Map present-pathogen indices into the full antigen vector
  array[nPp] int pres_idx;
  {
    int k = 1;
    for (p in 1:nP) {
      if (pres[p] == 1) {
        pres_idx[k] = p;
        k += 1;
      }
    }
  }
}

parameters {
  vector<lower=0, upper=1>[nPp] sero;   // infection prevalence per present pathogen
  vector<lower=0>[nP] sd0;              // per-antigen SD for negatives
  vector<lower=0>[nPp] sd1;            // per-pathogen SD for positives (v2 change)
  vector[nP] mu0;                       // mean for negatives (all antigens)
  vector<lower=0>[nPp] mu1;            // mean increment for present pathogens
  vector<lower=0>[nPp * (nP - 1)] phi; // cross-reactivity multipliers
  real<lower=0, upper=1> rho00;         // correlation among negatives
}

transformed parameters {
  array[nC] vector[nP] mu;                     // Gaussian means per component
  array[nC] vector[nP] sigma;                  // Gaussian SDs per component
  array[nC] cholesky_factor_cov[nP] L;         // Cholesky factors for fast likelihood
  matrix[N, nC] pC;                            // log p(y_n, c) per individual & component
  vector[N] log_lik;                           // individual log-likelihoods
  matrix[nPp, nP] CR;                         // cross-reactivity matrix (present × all)
  vector[nC] log_theta;                        // log mixture weights

  {
    array[nC] matrix[nP, nP] covM;             // local covariance matrices
    int ix = 1;

    // ---- 1. Cross-reactivity matrix ----
    for (p in 1:nPp) {
      for (a in 1:nP) {
        if (a == p) {
          CR[p, a] = 0;
        } else {
          CR[p, a] = phi[ix];
          ix += 1;
        }
      }
    }

    // ---- 2. Mixture weights (no strata) ----
    {
      vector[nP] log_sero_full   = rep_vector(0, nP);
      vector[nP] log1m_sero_full = rep_vector(0, nP);
      for (k in 1:nPp) {
        log_sero_full[pres_idx[k]]   = log(sero[k]);
        log1m_sero_full[pres_idx[k]] = log1m(sero[k]);
      }
      log_theta = infM * log_sero_full + infM_inv * log1m_sero_full;
    }

    // ---- 3. Component means & SDs ----
    // Component 1: all negative
    sigma[1] = sd0;
    mu[1]    = mu0;

    for (c in 2:nC) {
      // Initialise from all-negative baseline
      sigma[c] = sd0;
      mu[c]    = mu0;

      // Override positives: each uses its OWN sd1[a_pos]
      for (j in 1:npos[c]) {
        int a_pos = wpos[c, j];
        if (a_pos > 0) {
          sigma[c, a_pos] = sd1[a_pos];         // <-- per-pathogen sd1
          mu[c, a_pos]    = mu0[a_pos] + mu1[a_pos];
        }
      }

      // Override negatives with cross-reactivity contributions
      for (j in 1:(nP - npos[c])) {
        int a_neg = wneg[c, j];
        if (a_neg > 0) {
          real m = mu0[a_neg];
          real v = sd0[a_neg] ^ 2;               // per-antigen baseline variance
          for (k in 1:npos[c]) {
            int a_pos = wpos[c, k];
            m += CR[a_pos, a_neg] * mu1[a_pos];
            v += (CR[a_pos, a_neg] * sd1[a_pos]) ^ 2;  // <-- per-pathogen sd1
          }
          mu[c, a_neg]    = m;
          sigma[c, a_neg] = sqrt(v);
        }
      }
    }

    // ---- 4. Covariance matrices & Cholesky decomposition ----
    for (c in 1:nC) {
      // Diagonal: variances
      for (p in 1:nP) covM[c, p, p] = sigma[c, p] ^ 2;

      // Off-diagonal: covariances
      for (p in 1:(nP - 1)) {
        for (p2 in (p + 1):nP) {
          real cv;

          if (infM[c, p] + infM[c, p2] == 2) {
            // Both positive: independent
            cv = 0;
          } else if (infM[c, p] + infM[c, p2] == 0) {
            // Both negative: baseline correlation + cross-reactivity
            cv = rho00 * sd0[p] * sd0[p2];
            for (k in 1:npos[c]) {
              int a_pos = wpos[c, k];
              cv += CR[a_pos, p] * CR[a_pos, p2] * sd1[a_pos] ^ 2;  // <-- per-pathogen
            }
          } else {
            // One positive, one negative
            if (infM[c, p] == 1) {
              cv = CR[p, p2] * sd1[p] ^ 2;       // <-- per-pathogen
            } else {
              cv = CR[p2, p] * sd1[p2] ^ 2;      // <-- per-pathogen
            }
          }

          covM[c, p, p2] = cv;
          covM[c, p2, p] = cv;
        }
      }

      // Cholesky factor for efficient multivariate normal evaluation
      L[c] = cholesky_decompose(covM[c]);
    }

    // ---- 5. Vectorised likelihood via Cholesky ----
    for (c in 1:nC) {
      real constant_part = -0.5 * (nP * log_2pi + 2 * sum(log(diagonal(L[c]))));
      matrix[nP, N] z = mdivide_left_tri_low(L[c], Y - rep_matrix(mu[c], N));

      for (n in 1:N) {
        pC[n, c] = log_theta[c] + constant_part - 0.5 * dot_self(z[ : , n]);
      }
    }

    for (n in 1:N) log_lik[n] = log_sum_exp(pC[n,  : ]);
  }
}

model {
  // Priors (configurable via data)
  // NB: mu0/sd0 loc are vector[nP]; mu1/sd1 loc are vector[nPp]; scales scalar.
  //     `vector ~ normal(vector, real)` applies the prior element-wise.
  sero  ~ beta(sero_alpha, sero_beta);
  mu0   ~ normal(mu0_prior_loc, mu0_prior_scale);
  mu1   ~ normal(mu1_prior_loc, mu1_prior_scale);
  sd0   ~ normal(sd0_prior_loc, sd0_prior_scale);
  sd1   ~ normal(sd1_prior_loc, sd1_prior_scale);
  phi   ~ exponential(phi_rate);
  rho00 ~ beta(rho_alpha, rho_beta);

  target += sum(log_lik);
}

generated quantities {
  real sumloglik = sum(log_lik);

  // Posterior probability of each component for each individual
  array[N] simplex[nC] post_prob;
  for (n in 1:N) {
    post_prob[n] = softmax(to_vector(pC[n,  : ]));
  }
}
