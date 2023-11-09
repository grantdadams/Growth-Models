# This code simulates length-at-age data for individuals from a population where there are multiple sub-groups
# based on a von Bertalanfy growth curve

# The data are then fit to a hierarchical von bertalanfy growth curve with the following structure:

# Length.obs.i = Linf.i * (1-exp(-K.i * (age.obs.i - t0.i))) + eps.i
# Linf_i = exp(mu.linf + linf.group)
# K_i  = exp(mu.k + k.group)
# t0_i = mu.t0 + t0.group

# The group level parameters are assumed to be MVN with mean 0
# Eps.i is the normally distributed error term


####################################################################################
# Simulate length-at-age data ######################################################
####################################################################################
# Set Seed
library(MASS)
library(runjags)
library(ggplot2)
library(dplyr)
set.seed(1234)

## Specify data size
G=6
nsamples = 300
group <- sample(1:6, nsamples, replace = TRUE) # Group for individual X
N <- nsamples # Total number of samples


## Mu VBGM hyperparameters 
mu.Linf = 500
mu.k = 0.3 
mut.t0 = -0.5
mu.parms <- c(mu.Linf, mu.k, mut.t0)
sigma = 10 # Observation error


## Group level random effects
sigma.group = c(0.1, 0.05, 0.2)
rho = 0.3 # Correlation between group level parameters
cor.group.mat = matrix(rho, 3, 3)
diag(cor.group.mat) <- 1
cov.group.mat <- diag(sigma.group) %*% cor.group.mat %*% diag(sigma.group) # Get covariance


## Simulate parameters for groups----
# - Empty matrix and vectors to fill with parameters and data, respectively
group.param.mat <- group.re.mat <- matrix(NA,G,3,byrow = T)

# - Random effects
colnames(group.re.mat) <- c("log.Linf.group.re", "log.k.group.re", "t0.group.re")

# - On VBGF scale
colnames(group.param.mat) <- c("Linf.group", "k.group", "t0.group")


# - Simulate group level parameters
for(i in 1:G){
  # - Sim from mvnorm
  group.re.mat[i,] <- mvrnorm(1, rep(0,3), cov.group.mat) 
  
  # - Convert to parameter space
  group.param.mat[i,1:2] <- mu.parms[1:2] * exp(group.re.mat[i,1:2]) # Log to natural scale
  group.param.mat[i,3] <- mu.parms[3] + group.re.mat[i,3]
}


## Simulate length-at-age data ----
ages = seq(from=1,to=20, by = .05)
age = c()
length = c()
for(j in 1:N) {
  age[j] = sample(ages, 1) # Sample random age from age range
  length[j] = (group.param.mat[group[j],1] * (1 - exp(-group.param.mat[group[j],2]*(age[j]-group.param.mat[group[j],3])))) + rnorm(1,0,sigma) # Add normally distributed random error
  #FIXME: may want lognormal error
}


# Assign data to data frame
dat = data.frame(age = age, length = length, group = as.factor(group))
dat <- dat[which(dat$length > 0),] # Make sure all lengths are positive
dat <- dat %>% arrange(group, age, length)

# Plot the data
cols <- c("#86BBD8","#2F4858", "#F6AE2D", "#F26419", "#E86A92", "#57A773") # Colors for the VBGF lines (Females/Males)
ggplot(dat, aes(x = age, y = length, colour = group)) +
  geom_point(size = 2) + 
  scale_colour_manual(values=cols)


#################################################################################### 
# Stan VBGM ######################################################################## 
####################################################################################
## Stan setup ----
library(rstan)
rstan_options(threads_per_chain = 1)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)


# - Assign data to list
dat_list = list(
  Nobs = length(length),     # Number of obs
  G = G,                     # Number of groups
  Nages = ceiling(max(age)), # Max age (rounded) for predicting mean length-at-age by group
  length = length,           # Observed length
  age = age,                 # Age
  group = group,             # Group
  Zero = rep(0, 3)           # Vector of zero for mean of multivariate normal prior
)


## Fit model ----
# - VBGF with group level random effects
fit_stan <- stan(
  file = "R/5_hierarchical_growth_model_w_group_effects.stan",  # Stan program
  data = dat_list,             # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2000,          # number of warmup iterations per chain
  iter = 4000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  control = list(adapt_delta=0.99, stepsize=0.001, max_treedepth=18)
)
# - Some divergent transitions, may want to adjust


## Model diagnostics ----
# - Print and plot MCMC
# -- Look at chains
print(fit_stan, pars=c("mu_linf", "mu_k", "mu_t0", "linf_group", "k_group", "t0_group"), probs=c(.1,.5,.9)) # None of the betas are sig
traceplot(fit_stan, pars = c("mu_linf", "mu_k", "mu_t0", "sigma"), inc_warmup = FALSE, nrow = 2)
traceplot(fit_stan, pars = c("Lcorr_group", "sigma_group"), inc_warmup = FALSE, nrow = 2)

# -- Look at pairs
pairs(fit_stan, pars = c("mu_linf", "mu_k", "mu_t0"), las = 1)
pairs(fit_stan, pars = c("Lcorr_group", "sigma_group"), las = 1)

# - Sampler issues for all chains combined
sampler_params <- get_sampler_params(fit_stan, inc_warmup = TRUE)
summary(do.call(rbind, sampler_params), digits = 2)


## Plot model ----
# - Plot fitted model
plot(y = length , x = age, ylab = "Length", xlab = "Age", cex = 2, cex.lab = 1.25, 
     col = cols[group], pch = 16)
draws <- as.data.frame(fit_stan) # Get draws

# - Plot median of global curve
lines(1:max(age), apply(draws[,grepl("length_pred\\[1,",colnames(draws))], 2, median), col = 1, lty = 1, lwd = 4) # Global

# - Plot median of group curve
for(g in 1:G){
  # - Predicted
  lines(x = 1:max(ceiling(age)), 
        y = apply(draws[,grepl(paste0("length_pred\\[",g+1),colnames(draws))], 2, median), 
        col = cols[g], lty = 1, lwd = 2)
  
  # - True
  lines(x = 1:max(ceiling(age)), 
        y = (group.param.mat[g,1] * (1 - exp(-group.param.mat[g,2]*(1:max(ceiling(age))-group.param.mat[g,3])))), 
        col = cols[g], lty = 2, lwd = 2)
}

legend("bottomright", c("Global", paste0("Group ", 1:6)), col = c(1,cols), lty = 1, bty = "n", lwd = 2)


## Test group Parameters ----
# - Do the posteriers of the difference cross 0?
# - Only testing between group 1 and 2 here
par(mfrow = c(1,3))
hist(draws$`linf_group[1]`-draws$`linf_group[2]`, xlab = "Linf diff", main = NA)
hist(draws$`k_group[1]`-draws$`k_group[2]`, xlab = "K diff", main = "Group 1 and 2 difference")
hist(draws$`t0_group[1]`-draws$`t0_group[2]`, xlab = "t0 diff", main = NA)
par(mfrow = c(1,1))

