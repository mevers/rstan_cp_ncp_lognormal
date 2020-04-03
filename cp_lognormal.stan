data {
  int<lower=1> N;
  int<lower=1> J;
  vector<lower=0>[N] y;
  int<lower=1,upper=J> grp[N];
}

parameters {
  vector[J] theta;
  real<lower=0> sigma;

  // Hyperparameters
  real mu_theta;
  real<lower=0> sigma_theta;
}

transformed parameters {
  real log_sigma_theta = log(sigma_theta);
}

model {
  // Partial pooling
  theta ~ normal(mu_theta, sigma_theta);
  sigma ~ cauchy(0, 2.5);

  // Priors on the Hyperparameters
  mu_theta ~ normal(mean(log(y)), 5);
  sigma_theta ~ cauchy(0, 2.5);

  for (i in 1:N) {
    y[i] ~ lognormal(theta[grp[i]], sigma);
  }
}
