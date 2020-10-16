// Model code for von Bertalanffy growth function where growth parameters are conditionally autoregressive (CAR) random effects. Model 3; Eq. 13 from:
// Adams, G.D., Leaf, R.T., Ballenger, J.C., Arnott, S.A., Mcdonough, C.J., 2018. Spatial variability in the growth of Sheepshead (Archosargus probatocephalus) in the Southeast US : Implications for assessment and management. Fish. Res. 206, 35–43. doi:10.1016/j.fishres.2018.04.023


data {
  // ------------------------------------------------------------------------------------ //
  // 1. Assign data to stan objects                                                       //
  // ------------------------------------------------------------------------------------ //
  int<lower = 1>      n_i;                // Sample size
  int<lower = 1>      n_r;                // Number of regions
  int<lower = 1>      n_pred;             // Number of predictors
  vector[n_i]         log_length_i;       // Vector of log fork length of fish i
  vector[n_i]         age_i;              // Vector of age of fish i
  matrix[n_i, n_pred] design_mat;         // n_i * n_pred design matrix of linear predictors
  int<lower = 0>      r_i[n_i];           // Integer vector of region of fish i
  matrix<lower = 0,   upper = 1>[n_r, n_r] W; // Weights matrix where W_ij = 1 when regions i and j are neighbors
  matrix<lower = 0>   [n_r, n_r] D;       // Matrix where diagonal is number of neighbors of region i; diag(W * 1) (Eq. 13)
}
parameters {
  // ------------------------------------------------------------------------------------ //
  // 2. Specify model parameters                                                          //
  // ------------------------------------------------------------------------------------ //
  // 2.1. Regression coefficients for mean of VBGF parameters
  vector[n_pred] B_log_linf;
  vector[n_pred] B_log_k;
  vector[n_pred] B_t0;
  
  // 2.2. Random error (SD)
  real<lower=0> sigma;
  
  // 2.3. Hierarchical VBGF parameter error (SD)
  real<lower = 0> eta_linf;
  real<lower = 0> eta_k;
  real<lower = 0> eta_t0;
  
  // 2.4. Scaling parameters for non-centered distribution
  vector[n_r] alpha_linf;
  vector[n_r] alpha_k;
  vector[n_r] alpha_t0;
  
  // 2.5. Spatial-correlation factors
  real<lower = -1, upper = 1> phi_linf;
  real<lower = -1, upper = 1> phi_k;
  real<lower = -1, upper = 1> phi_t0;
}
transformed parameters{
  // ------------------------------------------------------------------------------------ //
  // 3. Specify derived model parameters to be saved                                      //
  // ------------------------------------------------------------------------------------ //
  // 3.1. Inverse of CAR model matrix
  matrix[n_r, n_r] L_linf;
  matrix[n_r, n_r] L_k;
  matrix[n_r, n_r] L_t0;
  
  // 3.2. VBGF parameter-specific random effects from CAR model
  vector[n_r] log_linf_re;
  vector[n_r] log_k_re;
  vector[n_r] t0_re;
  
  // 3.3. Precision of hierearchical error
  real tau_linf;
  real tau_k;
  real tau_t0;
  
  // 3.4. Get precision
  tau_linf  = 1/(eta_linf^2);
  tau_k     = 1/(eta_k^2);
  tau_t0    = 1/(eta_t0^2);

  // 3.5. Get inverse of CAR model matrix
  L_linf  = tau_linf * (D - phi_linf * W);
  L_k     = tau_k * (D - phi_k * W);
  L_t0    = tau_t0 * (D - phi_t0 * W);
  
  // 3.6. Get random effects
  log_linf_re = L_linf \ alpha_linf; // Equivalent to inverse(L) * alpha
  log_k_re    = L_k \ alpha_k;
  t0_re       = L_t0 \ alpha_t0;
}
model {
  // ------------------------------------------------------------------------------------ //
  // 4. Specify likelihood and priors                                                     //
  // ------------------------------------------------------------------------------------ //
  // 4.1. Temporary model objects to save parameter vectors
  vector[n_i] log_linf_i ;
  vector[n_i] k_i ;
  vector[n_i] t0_i ;
  
  // 4.2. Priors
  // 4.2.1. -- Regression coefficients for mean of VBGF parameters
  B_log_linf  ~ normal(0,10);  
  B_log_k     ~ normal(0,10);
  B_t0        ~ normal(0,10);
  
  // 4.2.2 -- Error components
  sigma       ~ cauchy(0, 5);
  eta_linf  ~ cauchy(0, 5);
  eta_k     ~ cauchy(0, 5);
  eta_t0    ~ cauchy(0, 5);
  
  // 4.2.3 -- Scaling factors for non-centered parameterization
  alpha_linf  ~ normal(0,1);
  alpha_k     ~ normal(0,1);
  alpha_t0    ~ normal(0,1);
  
  // 4.2.4. -- Spatial correlation parameters 
  phi_linf    ~ normal(0,10);
  phi_k       ~ normal(0,10);
  phi_t0      ~ normal(0,10);
  
  // 4.3. Model specification
  // 4.3.1. -- VBGF Parameters
  log_linf_i = design_mat * B_log_linf + log_linf_re[r_i];
  k_i        = exp(design_mat * B_log_k + log_k_re[r_i]);
  t0_i       = design_mat * B_t0 + t0_re[r_i]; 
  
  // 4.3.2. -- Model likelihood
  log_length_i ~ normal( log_linf_i + log1m_exp(-k_i .* (age_i - t0_i)), sigma);
}
