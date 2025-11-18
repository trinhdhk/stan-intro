// generated with brms 2.23.0
functions {
  /* helper function for asym_laplace_lpdf
   * Args:
   *   y: the response value
   *   quantile: quantile parameter in (0, 1)
   */
   real rho_quantile(real y, real quantile) {
     if (y < 0) {
       return y * (quantile - 1);
     } else {
       return y * quantile;
     }
   }
  /* asymmetric laplace log-PDF for a single response
   * Args:
   *   y: the response value
   *   mu: location parameter
   *   sigma: positive scale parameter
   *   quantile: quantile parameter in (0, 1)
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real asym_laplace_lpdf(real y, real mu, real sigma, real quantile) {
     return log(quantile * (1 - quantile)) -
            log(sigma) -
            rho_quantile((y - mu) / sigma, quantile);
   }
  /* asymmetric laplace log-CDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: location parameter
   *   sigma: positive scale parameter
   *   quantile: quantile parameter in (0, 1)
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real asym_laplace_lcdf(real y, real mu, real sigma, real quantile) {
     if (y < mu) {
       return log(quantile) + (1 - quantile) * (y - mu) / sigma;
     } else {
       return log1m((1 - quantile) * exp(-quantile * (y - mu) / sigma));
     }
   }
  /* asymmetric laplace log-CCDF for a single quantile
   * Args:
   *   y: a quantile
   *   mu: location parameter
   *   sigma: positive scale parameter
   *   quantile: quantile parameter in (0, 1)
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real asym_laplace_lccdf(real y, real mu, real sigma, real quantile) {
     if (y < mu) {
       return log1m(quantile * exp((1 - quantile) * (y - mu) / sigma));
     } else {
       return log1m(quantile) - quantile * (y - mu) / sigma;
     }
   }
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int<lower=1> NC_1;  // number of group-level correlations
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
}
transformed parameters {
  real quantile = 0.9;  // quantile parameter
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  // prior contributions to the log posterior
  real lprior = 0;
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
  lprior += student_t_lpdf(Intercept | 3, 5.1, 2.5);
  lprior += student_t_lpdf(sigma | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += student_t_lpdf(sd_1 | 3, 0, 2.5)
    - M_1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + Xc * b;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n];
    }
    for (n in 1:N) {
      target += asym_laplace_lpdf(Y[n] | mu[n], sigma, quantile);
    }
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
