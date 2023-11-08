# This script uses TMB to estimate length-at-age of a simulated population using a three parameter Von Bertalanffy Growth Function fit using maximum likelihood.

require(TMB)
# NOTE: You need RTools to make TMB function!!!!
# Follow the online instructions here:
# https://github.com/kaskr/adcomp

###### Simulate data
N_obs = 1000 # Number of observations
Min_age = 0 # Minimum age in the population
Max_age = 10 # Maximum age in the population

# Parameters
Obs_sd = 20 # Standard deviation
Linf = 350 # Asymptoptic length/mean length at maximum age (mm)
k = .3 # Growth rate (yr^-1)
t0 = -1.5 # Estimated age at length 0 (yr)

# Simulation
Age_i = runif(N_obs, Min_age, Max_age)
Obs_i = rnorm(N_obs, Linf * (1 - exp(-k * (Age_i - t0))), Obs_sd)

# Plot length-at-age
plot(Age_i,Obs_i , xlab = "Age (yr)", ylab = "Length (mm)")


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
int n = Age.size(); // get number of data points to loop over

// parameters:
PARAMETER(Linf); // asymptoptic length
PARAMETER(Kappa); // growth rate
PARAMETER(T_zero); // Age at length 0
PARAMETER(LogSigma); // log(residual SD)
// Fit sigma on a log scale to keep it > 0

// procedures: (transformed parameters)
Type sigma = exp(LogSigma);
vector LengthPred(n);

// Report transformed variables
ADREPORT(sigma);

Type nll = 0.0; // Initialize negative log-likelihood

for(int i = 0; i < n; i++){ // C++ starts loops at 0!
// get negative log likelihood (last argument is log = TRUE)
LengthPred(i) = Linf * (1.0-exp(-Kappa*(Age(i) - T_zero))) ;
}

nll = -sum(dnorm(Length,LengthPred , sigma,true)); // Calculate the negative log-likelihood
// NOTE: the negative log-likelihood can also be included in the for loop using the subtraction assignment '-=' where:
// nll -= dnorm(Length,LengthPred , sigma,true)
// This will subtract the each log-likelihood from the initial 0 value
return nll;
}"
# Write model to C++
write(tmb_model, file = "vbgf.cpp")

# Load model template
compile("vbgf.cpp") # Compile a C++ templated into a shared object file
dyn.load(dynlib("vbgf"))

# Run the model
obj = MakeADFun(
  data = list(Length = Obs_i, Age = Age_i), 
  parameters = list(Linf = 400, Kappa = .3,T_zero = 0 ,LogSigma = 0), DLL = "vbgf") # Prepare model object, data, and initial parameter estimates for fitting
opt1 = nlminb(obj$par, obj$fn, obj$gr) # Fit model using nlminb
summ = summary(sdreport(obj)) # Get summary of parameter estimates
summ

# Plot the fit
curve(summ[1,1] * (1 - exp(-summ[2,1] * (x - summ[3,1]))), from = 0 , to = Max_age, col = 2, add=T)