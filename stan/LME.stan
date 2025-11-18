// Demo stan code - LME model
// trinhdhk
functions {
    /* Row sum
    Row-wise summation of a matrix
    Args:
    * x: matrix
    Return: vector
    */
    vector row_sums(matrix x){
        return x * rep_vector(1, dims(x)[2]);
    }
}
data {
    // metadata
    int<lower=0> N; // sample size
    int<lower=0> M; // feature dimension
    int<lower=0> K; // random effect dimension
    int<lower=1> S; // number of group level

    // design matrix
    vector[N] Y;    // outcome
    matrix[N, M] X; // model matrix, including intercept
    matrix[N, K] Z; // model matrix for random effect, including intercept
    vector<lower=1, upper=S>[N] J_S;  // group level indicator
}

parameters {
    vector[M] beta;      // coef
    real<lower=0> sigma; // dispersion term
    vector<lower=0> sd;  // random effect dispersion term
    cholesky_factor_corr[K] L; // Lower triangular Cholesky decomposed correlation matrix
    matrix[S, K] zb;         // unstructured random effect
}

transformed parameters {
   matrix[S, K] b = transpose(diag_pre_multiply(sd, L) * zb); // random effect
   real lprior = 0.0; // if you want to investigate prior impact
   lprior += student_t_lpdf(beta | 3, 0, 2.5);
   lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5); //same as sigma ~ student_t(3,0,2.5)T[0, ];
   lprior += student_t_lpdf(sd | 3, 0, 2.5) - K * student_t_lccdf(0 | 3, 0, 2.5); //same as sd ~ student_t(3,0,2.5)T[0, ];
   lprior += lkj_corr_cholesky_lpdf(L | 1); // eta = 1 equivalent to Uniform(0,1);
}

model {
    vector[N] mu = X' * beta + row_sums(Z .* b[J_S,:]); // quotation mark for transpose(X) and * is matmul.
    target += lprior;       
    target += normal_lpdf(Y | mu, sigma);     // Y ~ normal(mu, sigma);
    target += std_normal_lpdf(to_vector(zb)); // to_vector(zb) ~ std_normal();
}

generated quantities {
   vector[N] log_lik;
   vector[N] mu = X' * beta + row_sums(Z .* b[J_S,:]);
   for (n in 1:N){
        log_lik[n] = normal_lpdf(Y[n] | mu[n], sigma);
   }
} 
