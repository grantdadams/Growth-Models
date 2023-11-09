//  VBGF with group level random effects
data {
  // Length-at-age data
  int<lower=0> Nobs;                    // number of observations 
  int<lower=0> G;                       // number of Groups 
  int<lower=0> Nages;                   // number of ages 
  real length[Nobs];                    // length
  real age[Nobs];                       // length
  int<lower=1, upper=G> group[Nobs];    // Sample group
  
  vector[3] Zero;                       // Vector of zero for mean of multivariate normal
}
parameters {
  // VBGF Params ----
  real mu_linf;                         // asymptotic length
  real mu_k;                            // growth coef
  real mu_t0;                           // age at length 0
  matrix[G, 3] eta_group;               // Group level deviation
  
  // Variance ----
  cholesky_factor_corr[3] Lcorr_group;  // Prior correlation for group-level variation
  vector<lower=0>[3] sigma_group;       // Prior scale for group-level variation
  
  real<lower=0> sigma;                  // observation error
}
transformed parameters {
  // Predicted length
  vector[Nobs] length_hat;
  
  // Group level parameters
  vector[G] linf_group = mu_linf * exp(eta_group[,1]);
  vector[G] k_group = mu_k * exp(eta_group[,2]);
  vector[G] t0_group = mu_t0 + eta_group[,3];
  
  // Predicted length
  for(i in 1:Nobs){
    length_hat[i] = linf_group[group[i]] * (1-exp(-k_group[group[i]] * (age[i] - t0_group[group[i]])));
  }
}
model {
  // Priors
  // - Global parameters (will likely adjust these based on species)
  mu_linf ~ lognormal(log(500), 0.2);
  mu_k ~ lognormal(log(0.3), 0.2);
  mu_t0 ~ normal(0, 0.5);
  
  // - Group level variation priors
  sigma_group ~ cauchy(0, 0.5);
  Lcorr_group ~ lkj_corr_cholesky(1);
  
  for(i in 1:G){
    eta_group[i,] ~ multi_normal_cholesky(Zero, diag_pre_multiply(sigma_group, Lcorr_group));
  }
  
  // Likelihood
  length ~ normal(length_hat, sigma);
} 
generated quantities{
  // Predicted length
  matrix[1+G, Nages] length_pred;
  
  // Global
  for(i in 1:Nages){
    length_pred[1,i] = mu_linf * (1 - exp(-mu_k * (i - mu_t0))); // Global
  }
  
  // Loop through groups
  for(g in 1:G){
    // Loop through ages
    for(i in 1:Nages){
      length_pred[g+1,i] = mu_linf * exp(eta_group[g,1]) * (1 - exp(-mu_k * exp(eta_group[g,2]) * (i - mu_t0 - eta_group[g,3]))); 
    }
  }
}
