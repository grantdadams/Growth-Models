## SETUP ----
# Install packages (for rethinking)
install.packages(c("coda","mvtnorm","devtools","loo","dagitty","shape", "cmdstanr", "invgamma", "ggdist"))
# - May need this for cmdstanr
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
devtools::install_github("rmcelreath/rethinking")

# Load packages and data
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]

ndraws <- 1000
nobs <- length(d2$weight)


## PRIOR MODEL 1 & 2 ----
# - Set priors for default models (page 93)
# -- Model 1 = intercept only
# -- Model 2 = Effect of weight
alpha_prior <- rnorm(ndraws, mean = 178, sd = 100) # Intercept prior
beta_prior <- rnorm(ndraws, mean = 0, sd = 10) # Regression prior for weight
sd_prior <- runif(ndraws, 0, 50)

# - Prior predictive distribution of height
prior_predictive_height_mod1 <- matrix(NA, nrow = nobs, ncol = ndraws)
prior_predictive_height_mod2 <- matrix(NA, nrow = nobs, ncol = ndraws) 

# -- Can be vectorized for speed, but this illustrates the process better
for(i in 1:ndraws){ # Loop through prior draws
  for(j in 1:nobs){ # Loop through data
    prior_predictive_height_mod1[j,i] <- rnorm(1, mean = alpha_prior[i], sd = sd_prior[i])
    prior_predictive_height_mod2[j,i] <- rnorm(1, mean = alpha_prior[i] + beta_prior[i] * d2$weight[j], sd = sd_prior[i])
  }
}


## PRIOR MODEL 3 ----
# Model with regularized horseshoe prior on beta for effect of weight
# - Horseshoe prior is a form of regularization that pulls regression coefficients towards zero
# -- See https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-11/issue-2/Sparsity-information-and-regularization-in-the-horseshoe-and-other-shrinkage/10.1214/17-EJS1337SI.full
# -- Appendix C has the code
prior_predictive_height_mod3 <- matrix(NA, nrow = nobs, ncol = ndraws)

# - Update beta prior to horseshoe prior
library(invgamma)
library(ggdist)

p0 <- 2 # prior guess for the number of relevant variables
D <- 4
tau0 <- p0 /( D - p0 ) / sqrt (nobs) # rstanarm will scale this by sigma automatically

# - Prior specifications
scale_global <- tau0 # scale for the half -t prior for tau
global_df <- 1 # degrees of freedom for the half -t prior for tau
local_df  <- 1 # degrees of freedom for the half - t priors for lambdas. Setting to 1 is horshoe
slab_scale <- 2.5 # slab scale for the regularized horseshoe
slab_df <- 4 # slab degree

# - Priors
z <-  rnorm(ndraws, mean = 0, sd = 1)
lambda <- rstudent(ndraws, local_df, 0, 1)
tau <- rstudent(ndraws, global_df, 0, scale_global * sd_prior)
caux <- rinvgamma(ndraws, 0.5* slab_df , 0.5* slab_df )

# -- Alternative parameterization
# aux1_local <-  rnorm(ndraws, mean = 0, sd = 1)
# aux2_local <- rinvgamma(ndraws, 0.5* local_df , 0.5* local_df)
# aux1_global <-  rnorm(ndraws, mean = 0, sd = 1)
# aux2_global <- rinvgamma(ndraws, 0.5* global_df , 0.5* global_df)
# 
# lambda = aux1_local * sqrt( aux2_local )
# tau =  aux1_global * sqrt ( aux2_global ) * scale_global * sd_prior

# - Derived parameters
c = slab_scale * sqrt( caux )
lambda_tilde = sqrt ( c ^2 *( lambda^2 ) / (c ^2 + tau^2*( lambda^2 )) )
beta_horshoe_prior <- z * lambda_tilde*tau

for(i in 1:ndraws){
  for(j in 1:nobs){
    prior_predictive_height_mod3[j,i] <- rnorm(1, mean = alpha_prior[i] + beta_horshoe_prior[i] * d2$weight[j], sd = sd_prior[i])
  }
}


## PLOTS ----
# Plot distribution of betas
par(mfrow = c(1,2))
hist(beta_prior, xlab = "Beta (cm)", main = "Beta prior")
hist(beta_horshoe_prior, xlab = "Beta (cm)", main = "Regularized horseshoe beta prior")

# Plot marginal distribution of heights and prior predictive heights
par(mfrow = c(2,2))
hist(d2$height, xlab = "Observed height (cm)", main = "Observed")
hist(prior_predictive_height_mod1, xlab = "Prior predictive height (cm)", main = "Model 1: intercept only")
hist(prior_predictive_height_mod2, xlab = "Prior predictive height (cm)", main = "Model 2: +weight")
hist(prior_predictive_height_mod3, xlab = "Prior predictive height (cm)", main = "Model 3: +weight with horseshoe prior")
