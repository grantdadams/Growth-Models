// The input data
data {
  int<lower=0> nages;
  int<lower=0> nyrs;
  matrix[nyrs, nages] wt;
  
  # Survey data
  matrix[nsrvs, nyrs] srv_obs;
  real srv_comp_obs[nsrvs, nyrs, nages];
  
  # Fishery data
  vector[nyrs] catch_obs;
  matrix[nyrs, nages] catch_comp_obs;
}
// The parameters accepted by the model.
parameters {
  // Recruitment
  real<lower=0> log_mean_R;
  vector[nyrs] R_dev;
  real<lower=0> R_sigma;
  
  // Survey observation model
  real<lower=0> srv_slope[nsrvs];   // Logistic selectivity slope 
  real srv_asymp[nsrvs];            // Logistic selectivity asymptote 
  real srv_q[nsrvs];                // Survey catchability
  
  // Fishery observation model
  real<lower=0> fsh_slope;          // Logistic selectivity slope 
  real fsh_asymp;                   // Logistic selectivity asymptote 
  real F[nyrs];                     // Fishing mortality
}
// The age-structured model to be estimated
transformed parameters{
  matrix[nyrs, nages] n_at_age;     // Numbers-at-age
  matrix[nyrs, nages] f_at_age;     // Numbers-at-age
  
  vector[nyrs] catch_hat;                  // Predicted catch
  matrix[nyrs, nages] catch_comp_hat;      // Predicted fishery catch-at-age
  matrix[nsrvs, nyrs] srv_hat;             // Predicted sruvey
  real srv_comp_hat[nsrvs, nyrs, nages]; // Predicted survey catch-at-age
  
  // Recruitment
  n_at_age[,1] = exp(log_mean_R + R_dev);
  
  // Initial abundance (assuming starting at unfished equilibrium)
  for(age in 2:nages){
    n_at_age[1,age] = n_at_age[1, age-1]*exp(-M); // Could add initial deviates here
  }
  n_at_age[1,nages] /= (1.0 - exp(-M));
  
  // Calculate fishing mortality-at-age
  for(yr in 1:nyrs){
    for(age in 1:nages){
      f_at_age[yr, age] = 1/(1+exp(-fsh_slope * (age - fsh_asymp))) * F[yr];
    }
  }
  
  // Numbers at age (years > 1 and ages > 1)
  for(yr in 2:nyrs){
    for(age in 2:nages){
      n_at_age[yr,age] = n_at_age[yr-1,age-1] * exp(-M - F_at_age[yr-1, age-1]);
    }
  }
  
  // Survey data
  for(srv in 1:nsrvs){
    for(yr in 1:nyrs){
      srv_hat[srv, yr] = 0;
      for(age in 1:nages){
        srv_comp_hat[srv, yr, age] = srv_q[srv]/(1+exp(-srv_slope[srv] * (age - srv_asymp[srv])));
        srv_hat[srv, yr] += srv_comp_hat[srv, yr, age] * wt[yr, age];
      }
      
      // - Normalize
      srv_comp_hat[srv, yr, age] /= sum(srv_comp_hat[srv, yr,]);
    }
  }
}
// Model likelihoods
model{
  // Recruitment deviates
  R_dev ~ lognormal(log_mean_R, R_sigma);
  
  // Survey likelihoods
  
  // Fishery likelihoods
  y ~ normal(mu, sigma);
}

