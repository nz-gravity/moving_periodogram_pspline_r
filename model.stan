data {
  int<lower=1> N;
  int<lower=1> Kt;
  int<lower=1> Kf;
  int<lower=1> J;
  matrix[N, Kt] Bt;
  matrix[J, Kf] Bf;
  array[N] int<lower=1, upper=J> rung;
  vector<lower=0>[N] y;
  vector<lower=0>[Kt*Kf] lt;
  vector<lower=0>[Kt*Kf] lf;
  array[Kt*Kf] int<lower=0, upper=1> is_null;
  real<lower=0> null_precision;
  real<lower=0> ridge;
  real<lower=0> alpha_phi;
  real<lower=0> beta_phi;
}
parameters {
  vector[Kt*Kf] c;
  vector<lower=0>[2] phi;
}
model {
  matrix[Kt, J] cf = to_matrix(c, Kt, Kf) * Bf';
  phi ~ gamma(alpha_phi, beta_phi);
  for (k in 1:(Kt*Kf)) {
    real precision = is_null[k] ? null_precision : phi[1]*lt[k] + phi[2]*lf[k] + ridge;
    c[k] ~ normal(0, inv_sqrt(precision));
  }
  {
    vector[N] eta = rows_dot_product(Bt, cf[, rung]');
    target += -sum(eta) - dot_product(y, exp(-eta));
  }
}
