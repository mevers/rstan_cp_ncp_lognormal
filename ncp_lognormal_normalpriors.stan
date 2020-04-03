data {
  int<lower=1> N;
  int<lower=1> J;
  vector<lower=0>[N] y;
  int<lower=1,upper=J> grp[N];
}

parameters {
  vector[J] theta_raw;
  real<lower=0> sigma;

  // Hyperparameters
  real mu_theta;
  real<lower=0> sigma_theta;
}

transformed parameters {
  vector[J] theta;
  real log_sigma_theta = log(sigma_theta);

  // Non-centred parametrisation
  // This is the same as theta ~ normal(mu_d, sigma_d)
  for (j in 1:J) {
    theta = mu_theta + sigma_theta * theta_raw;
  }
}

model {
  // Prior on non-centred theta and sigma
  theta_raw ~ std_normal();
  sigma ~ normal(0, 5);

  // Priors on the Hyperparameters
  mu_theta ~ normal(mean(log(y)), 5);
  sigma_theta ~ normal(0, 5);

  for (i in 1:N) {
    y[i] ~ lognormal(theta[grp[i]], sigma);
  }
}
