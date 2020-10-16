# This code creates a 3-parameter von Bertalanffy growth function assuming a left and right truncated lognormal distribution in template model builder.

# Load packages
library(TMB)

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