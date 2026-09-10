functions {
  real chirp_loglik(real A, real cycles_f0, real cycles_fdot, vector c,
                    matrix B, vector t, array[] int index, int L,
                    matrix wr, matrix wi, vector dre, vector dimag) {
    int N = rows(B);
    vector[num_elements(t)] wave = A*cos(2*pi()*(cycles_f0*t + 0.5*cycles_fdot*square(t)));
    matrix[N,L] windows = to_matrix(wave[index],N,L);
    vector[N] rr = dre - rows_dot_product(windows,wr);
    vector[N] ri = dimag - rows_dot_product(windows,wi);
    vector[N] eta = B*c;
    return -sum(eta) - dot_product(square(rr)+square(ri),exp(-eta));
  }
}
data {
  int<lower=1> T;
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1> L;
  vector[T] t;  // zero-based time divided by T; dt=1 second
  array[N*L] int<lower=1,upper=T> index;
  matrix[N,L] wr;
  matrix[N,L] wi;
  vector[N] dre;
  vector[N] dimag;
  matrix[N,K] B;
  vector<lower=0>[K] lt;
  vector<lower=0>[K] lf;
  real f0_lower;
  real f0_upper;
  real fdot_lower;
  real fdot_upper;
}
parameters {
  real A;
  real<lower=f0_lower,upper=f0_upper> cycles_f0;
  real<lower=fdot_lower,upper=fdot_upper> cycles_fdot;
  vector[K] c;
  vector<lower=0>[2] phi;
}
model {
  A ~ normal(0,5);
  phi ~ gamma(2,1);
  // Uniform signal priors inside the specified cycle-coordinate bounds.
  for (k in 1:K) {
    real q = (lt[k]+lf[k] == 0) ? 0.01 : phi[1]*lt[k]+phi[2]*lf[k]+1e-6;
    c[k] ~ normal(0,inv_sqrt(q));
  }
  target += chirp_loglik(A,cycles_f0,cycles_fdot,c,B,t,index,L,wr,wi,dre,dimag);
}
generated quantities {
  real f0 = cycles_f0/T;          // Hz
  real fdot = cycles_fdot/(1.0*T*T); // Hz/s
  // Also makes an independent R-versus-Stan likelihood check straightforward.
  real log_lik = chirp_loglik(A,cycles_f0,cycles_fdot,c,B,t,index,L,wr,wi,dre,dimag);
}
