data {
  int<lower=1> N;                      // number of observations
  vector<lower=0>[N] Cyt;            // Cytosolic counts
  vector<lower=0>[N] Nuc;           // Nuclear counts
  vector<lower=0>[N] Who;          // Whole-cell counts
}

transformed data {
  real<lower=0> sigma_data;
  sigma_data = sd(Who);
}

parameters {
  real<lower=0, upper=1> beta;    // coefficient
  real<lower=0> sigma;            // scale of the t-distribution
  real<lower=1> nu;               // degrees of freedom of the t-distribution
}

model {
  // Likelihood
  // Student's t-distribution instead of normal for robustness
  vector[N] mu;
  mu = beta * Cyt + (1 - beta) * Nuc; // mean response
  Who ~ student_t(nu, mu, sigma);
  // Uninformative priors on all parameters
  beta ~ beta(2, 2);
  sigma ~ exponential(1/sigma_data);
  nu ~ gamma(2, 0.1);
}
