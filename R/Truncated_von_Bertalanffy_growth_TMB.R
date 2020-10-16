# This code creates a 3-parameter von Bertalanffy growth function assuming a left and right truncated normal distribution of log length in template model builder.

# NOTE: You need RTools to make TMB function!!!!
# Follow the online instructions here:
# https://github.com/kaskr/adcomp

###### Simulate data
N_obs = 1000 # Number of observations
Min_age = 0 # Minimum age in the population
Max_age = 10 # Maximum age in the population
min_size = 40 # Minimum size harvested by the fishery
max_size = 400 # Maximum size harvested by the fishery


# Parameters
Obs_sd = 20 # Standard deviation
Linf = 350 # Asymptoptic length/mean length at maximum age (mm)
k = .3 # Growth rate (yr^-1)
t0 = -1.5 # Estimated age at length 0 (yr)

# Simulation
Age_i = runif(N_obs, Min_age, Max_age)
Obs_i = rnorm(N_obs, Linf * (1 - exp(-k * (Age_i - t0))), Obs_sd)

# Truncate data
Age_i = Age_i[which(Obs_i >= min_size & Obs_i <= max_size)]
Obs_i = Obs_i[which(Obs_i >= min_size & Obs_i <= max_size)]

# Plot length-at-age
plot(Age_i,Obs_i , xlab = "Age (yr)", ylab = "Length (mm)")

#..........................................................................
# Write model
# Remember the language is C++ so comments are indicated by "//" and the end of an expression needs to be explicitly stated by typing ";"
tmb_model = "
// von Bertalanffy growth function
#include 

template
Type objective_function::operator() () {
// data:
DATA_VECTOR(Age);
DATA_VECTOR(Length);
DATA_SCALAR(log(min_size)); // Minimum size for truncation
DATA_SCALAR(log(max_size)); // Maximum size for truncation
int n = Age.size(); // get number of data points to loop over

// VBGF parameters:
PARAMETER(log_Linf); // asymptoptic length
PARAMETER(log_Kappa); // growth rate
PARAMETER(log_T_zero_plus_10); // Age at length 0
PARAMETER(log_Sigma); // log(residual SD)
// Fit sigma on a log scale to keep it > 0

// procedures: (transformed parameters)
Type Linf = exp(log_Linf);
Type Kappa = exp(log_Kappa);
Type T_zero = exp(log_T_zero_plus_10) - 10;
Type sigma = exp(log_Sigma);

// Report transformed variables
ADREPORT(sigma);
ADREPORT(Linf);
ADREPORT(Kappa);
ADREPORT(T_zero);

// Vectors of data and predictors
vector LengthPred(n);
vector CumulativeNorm(n);
vector TruncatedNorm(n);

Type nll = 0.0; // Initialize negative log-likelihood

for(int i = 0; i < n; i++){ // C++ starts loops at 0!

// Extended VBGF Parameters
Type KappaTemp = exp(log_Kappa);
Type LinfTemp = exp(log_Linf);

// VBGF
LengthPred(i) = log(LinfTemp * (1.0-exp(-KappaTemp*(Age(i) - T_zero))));

// Calculate the CDF for the right truncation minus the CDF for the left truncation
CumulativeNorm(i) = pnorm(max_size,LengthPred(i),sigma)  - pnorm(min_size,LengthPred(i),sigma);

// Calculate the truncated nromal distribution
TruncatedNorm(i) = (dnorm(Length(i),LengthPred(i) , sigma, false)) / CumulativeNorm(i);

// Calculate the negative log-likelihood
nll -= log(TruncatedNorm(i)); 
}
return nll;
}"
# Write model to C++
write(tmb_model, file = "vbgf_trunc.cpp")

# Load model template
compile("vbgf_trunc.cpp") # Compile a C++ templated into a shared object file
dyn.load(dynlib("vbgf_trunc"))

# Run the model
obj = MakeADFun(
  data = list(Length = Obs_i, Age = Age_i, min_size = min_size, max_size = max_size), 
  parameters = list(Linf = 400, Kappa = .3,T_zero = 0 ,LogSigma = 0), DLL = "vbgf_trunc") # Prepare model object, data, and initial parameter estimates for fitting
opt1 = nlminb(obj$par, obj$fn, obj$gr) # Fit model using nlminb
summ = summary(sdreport(obj)) # Get summary of parameter estimates
summ

# Plot the fit
curve(summ[1,1] * (1 - exp(-summ[2,1] * (x - summ[3,1]))), from = 0 , to = Max_age, col = 2, add=T)
