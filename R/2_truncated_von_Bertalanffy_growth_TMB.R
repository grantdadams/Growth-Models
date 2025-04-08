# This code creates a 3-parameter von Bertalanffy growth function assuming a left and right truncated normal distribution of log length in template model builder.

# NOTE: You need RTools to make TMB function!!!!
# Follow the online instructions here:
# https://github.com/kaskr/adcomp

###### Simulate data
library(RTMB)
N_obs = 1000 # Number of observations
Min_age = 0 # Minimum age in the population
Max_age = 10 # Maximum age in the population
min_size = 40 # Minimum size harvested by the fishery
max_size = 300 # Maximum size harvested by the fishery


# Parameters
Obs_sd = 20 # Standard deviation
Linf = 350 # Asymptoptic length/mean length at maximum age (mm)
k = .3 # Growth rate (yr^-1)
t0 = -1.5 # Estimated age at length 0 (yr)

# Simulation
Age_i = runif(N_obs, Min_age, Max_age)
Length_i = rnorm(N_obs, Linf * (1 - exp(-k * (Age_i - t0))), Obs_sd)

# Truncate data
Obs_age_i = Age_i[which(Length_i >= min_size & Length_i <= max_size)]
Obs_length_i = Length_i[which(Length_i >= min_size & Length_i <= max_size)]

# Plot length-at-age
plot(Age_i, Length_i , xlab = "Age (yr)", ylab = "Length (mm)", pch = 16)
points(Obs_age_i, Obs_length_i , col = 2, pch = 16)

#..........................................................................
# Write model
# Remember the language is C++ so comments are indicated by "//" and the end of an expression needs to be explicitly stated by typing ";"
vbgf_trunc = function(parList){
  
  # Data
  Age = dataList$Age
  Length = dataList$Length
  min_size = dataList$min_size
  max_size = dataList$max_size
  
  # VBGF parameters:
  Linf = exp(parList$log_Linf)    # Asymptoptic length
  Kappa = exp(parList$log_Kappa)  # growth rate
  T_zero = parList$T_zero         # Age at length 0
  Sigma = exp(parList$log_Sigma); # residual SD. Fit sigma on a log scale to keep it > 0
  
  # Model
  LengthPred = Linf * (1.0-exp(-Kappa*(Age - T_zero)))
  
  # Calculate the CDF for the right truncation minus the CDF for the left truncation
  CumulativeNorm = pnorm(max_size, LengthPred, Sigma) - pnorm(min_size, LengthPred, Sigma)
  
  # Calculate the truncated nromal distribution
  TruncatedNorm = (dnorm(Length, LengthPred , Sigma, FALSE)) / CumulativeNorm
  
  # Calculate the negative log-likelihood
  nll = -sum(log(TruncatedNorm))
  
  # Report transformed variables
  REPORT(LengthPred)
  REPORT(CumulativeNorm)
  REPORT(TruncatedNorm)
  REPORT(Sigma)
  REPORT(Linf)
  REPORT(Kappa)
  REPORT(T_zero)
  
  # Get SD of transformed variables
  ADREPORT(Sigma)
  ADREPORT(Linf)
  ADREPORT(Kappa)
  ADREPORT(T_zero)
  
  return (nll)
}

# Run the model
dataList <- list(Length = Obs_length_i, Age = Obs_age_i, min_size = min_size, max_size = max_size)
parList <- list(log_Linf = log(400), log_Kappa = log(.3), T_zero = 0 , log_Sigma = 3)

obj = RTMB::MakeADFun(
  vbgf_trunc,
  parameters = parList) # Prepare model object, data, and initial parameter estimates for fitting
opt1 = nlminb(obj$par, obj$fn, obj$gr) # Fit model using nlminb
summ = summary(sdreport(obj)) # Get summary of parameter estimates
summ

# Plot the fit
curve(summ[6,1] * (1 - exp(-summ[7,1] * (x - summ[8,1]))), from = 0 , to = Max_age, col = 3, lwd = 3, add=T)
