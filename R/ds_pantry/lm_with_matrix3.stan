data {
  int<lower=0> N; // Número de observaciones
  int<lower=0> K; // Número de predictores
  
  vector[N] y; // Outcome

  // row_vector[m] L; // Lower bounds
  // row_vector[m] U; // Upper bounds
  
  matrix[N, K] x; // Matriz de predictores
  
  // int intcpt; // Modelo con (= 1) o sin (= 0) intercept
}

transformed data {
  matrix[N, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  
  // thin and scale the QR decomposition
  Q_ast = qr_Q(x)[, 1:K] * sqrt(N - 1);
  R_ast = qr_R(x)[1:K, ] / sqrt(N - 1);
  R_ast_inverse = inverse(R_ast);
}


parameters {
  // real alpha_raw;           // intercept
  vector[K] theta;      // coefficients on Q_ast
  real<lower=0> sigma;  // error scale
}

transformed parameters {
  // real alpha;
  
  // alpha = alpha_raw * intcpt;
  
}


model {
  
  // LIKELIHOOD
  // y ~ normal(Q_ast * theta + alpha, sigma);  
  y ~ normal(Q_ast * theta, sigma); 
}

generated quantities {
  vector[K] beta;
  beta = R_ast_inverse * theta; // coefficients on x
}

