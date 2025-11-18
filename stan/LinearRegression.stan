// Demo stan code - Linear model
// trinhdhk
data {
    // metadata
    int<lower=0> N; // sample size
    int<lower=0> M; // feature dimension

    // design matrix
    vector[N] Y;    // outcome
    matrix[N, M] X; // model matrix, including intercept
}

parameters {
    vector[M] beta;      // coef
    real<lower=0> sigma; // dispersion term
}

transformed parameters {
   real lprior = 0.0; // if you want to investigate prior impact
   lprior += student_t_lpdf(beta | 3, 0, 2.5);
   lprior += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5); //same as sigma ~ student_t(3,0,2.5)T[0, ];
}

model {
    vector[N] mu = X' * beta; // quotation mark for transpose(X) and * is matmul.
    target += lprior;       
    target += normal_lpdf(Y | mu, sigma); 
}

generated quantities {
   vector[N] log_lik;
   vector[N] mu = X' * beta; 
   for (n in 1:N){
        log_lik[n] = normal_lpdf(Y[n] | mu[n], sigma);
   }
} 
