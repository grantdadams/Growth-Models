# This code simulates length-at-age data for individuals from a population where there are multiple sub-groups and there are sex and location (two rivers) specific differences in growth. Multiple length-at-age samples are "simulated/sampled" from each individual, representing capture-recapture data.

# The data are then fit to a hierarchical von bertalanfy growth curve with the following structure:

# Length.obs.i = Linf.i * (1-exp(-K.i * (age.obs.i - t0.i))) + eps.i
# Linf_i = exp(mu.linf + linf.group + linf.i + beta.linf.river * river.i + beta.linf.sex * sex.i)
# K_i  = exp(mu.k + k.group + k.i + beta.k.river * river.i + beta.k.sex * sex.i)
# t0_i = mu.t0 + t0.group + t0.i + beta.t0.river * river.i + beta.t0.sex * sex.i

# Where the "betas" are river and sex effects depending on the river of origin and sex of individual "i" 
# The group and and individual i level parameters are assumed to be MVN with mean 0 and conjugate prior variance
# Eps.i is the normally distributed error term
# NOTE: given the non-linearity and high dimension of the model, estimation is probably better in Stan


####################################################################################
# Simulate length-at-age data ######################################################
####################################################################################
# Set Seed
library(MASS)
library(runjags)
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
sigma.group = c(0.01, 0.05, 0.2)
rho = 0.3 # Correlation between group level parameters
cor.group.mat = matrix(rho, 3, 3)
diag(cor.group.mat) <- 1
cov.group.mat <- diag(sigma.group) %*% cor.group.mat %*% diag(sigma.group) # Get covariance


## Simulate parameters for groupsv----
# - Empty matrix and vectors to fill with parameters and data, respectively
group.param.mat <- group.re.mat <- matrix(NA,G,3,byrow = T)

# - Random effects
colnames(group.re.mat) <- c("log.Linf.group.re", "log.k.group.re", "t0.group.re")

# - On VBGF scalge
colnames(group.param.mat) <- c("Linf.group", "k.group", "t0.group")


# - Simulate group level parameters
for(i in 1:G){
  group.re.mat[i,] <- mvrnorm(1, rep(0,3), cov.group.mat)
  group.param.mat[i,1:2] <- mu.parms[1:2] * exp(group.re.mat[i,1:2]) # Log to natural scale
  group.param.mat[i,3] <- mu.parms[3] + group.re.mat[i,3]
}


## Simulate length-at-age data
ages = seq(from=1,to=20, by = .05)
age = c()
length = c()
for(j in 1:N) {
  age[j] = sample(ages, 1) # Sample random age from age range
  length[j] = (ind.param.mat[sample.id[j],1] * (1 - exp(-ind.param.mat[sample.id[j],2]*(age[j]-ind.param.mat[sample.id[j],3])))) + rnorm(1,0,sigma)
}


# Assign data to data frame
dat = data.frame(age = age, length = length, sample.group = sample.group, sample.id = sample.id, sample.sex = sample.sex, sample.river = sample.river)
dat <- dat[which(dat$length > 0),]
dat <- dat %>% arrange(sample.group, sample.id, age)

# Plot the data
for (i in 1:length(unique(dat$sample.group))){
  sub = dat[which(dat$sample.group==i),]
  if ( i == 1){
    plot(sub$age, sub$length, pch = sub$sample.sex+1 ,col = i, xlab = "Age", ylab = "Length", ylim = c(0, max(length, na.rm = TRUE)), xlim = range(ages))
  } else {
    points(sub$age, sub$length , col = i, pch = sub$sample.sex+1)
  }
}



#################################################################################### 
# JAGS VBGM ######################################################################## 
####################################################################################
require(runjags)
require(coda)

##### SPECIFY JAGS MODEL CODE #####
jags.mod =
  "model{

## Data likelihood ----
for(i in 1:N){
  length[i] ~ dnorm(y.hat[i], tau.y) # May want this to be lognormal
  y.hat[i] = Linf.i[sample.id[i]] * (1-exp(-k.i[sample.id[i]] * (age[i] - t0.i[sample.id[i]])))
}

# - Observation error
tau.y = pow(sigma, -2)
sigma ~ dunif(0,100)


## Group level random effects ----
for (j in 1:G){ # - Loop through groups
  group.re.pars[j,1:3] ~ dmnorm(mu.re, group.inv.cov.mat)  # The `0` may need to be a vector of 3 zeros?

  log.Linf.group.re[j] = group.re.pars[j,1]
  log.k.group.re[j] = group.re.pars[j,2] 
  t0.group.re[j] = group.re.pars[j,3]
  
  # - Get group level params on natural scale for reporting
  Linf.group[j] = exp(log.mu.Linf + log.Linf.group.re[j])
  k.group[j] = exp(log.mu.k + log.k.group.re[j])
  t0.group[j] = mu.t0 + t0.group.re[j]
}

# - Group level covariance matrix
group.inv.cov.mat ~ dwish( z.group.mat, 4 ) # Add `z.group.mat = diag(x=1,nrow=3)` to the data
group.cov.mat = inverse( group.inv.cov.mat ) # you can get out the variance and correlation coeficients using `cov2cor` in R or add code here


## Individual level random effects ----
for (k in 1:Nind){ # - Loop through individuals
  ind.re.pars[k,1:3] ~ dmnorm(mu.re, ind.inv.cov.mat) # The `0` may need to be a vector of 3 zeros?

  log.Linf.ind.re[k] = ind.re.pars[k,1]
  log.k.ind.re[k] = ind.re.pars[k,2] 
  t0.ind.re[k] = ind.re.pars[k,3]
  
  # - Get individual level params on natural scale 
  # NOTE: group is now of length Nind and will have to be adjusted in the data
  Linf.i[k] = exp(log.mu.Linf + log.Linf.group.re[group[k]] + log.Linf.ind.re[k] + log.linf.beta.sex * sex[k] + log.linf.beta.river * river[k]) # Exponent here to keep Linf positive
  k.i[k] = exp(log.mu.k + log.k.group.re[group[k]] + log.k.ind.re[k] + log.k.beta.sex * sex[k] + log.k.beta.river * river[k]) # Exponent here to keep K positive
  t0.i[k] = mu.t0 + t0.group.re[group[k]] + t0.ind.re[k] + t0.beta.sex * sex[k] + t0.beta.river * river[k]
}

# - Individual level covariance matrix
ind.inv.cov.mat ~ dwish( z.ind.mat, 4) # Add `z.ind.mat = diag(x=1,nrow=3)` to the data
ind.cov.mat = inverse( ind.inv.cov.mat ) # you can get out the variance and correlation coeficients using =`cov2cor` in R or add code here

# Sex and river effet priors
log.linf.beta.river ~ dnorm(0, 0.01)
log.k.beta.river ~ dnorm(0, 0.01)
t0.beta.river ~ dnorm(0, 0.01)

log.linf.beta.sex ~ dnorm(0, 0.01)
log.k.beta.sex ~ dnorm(0, 0.01)
t0.beta.sex ~ dnorm(0, 0.01)

## Population level parameters ----
# - These could probably be more informative given previous studies
log.mu.Linf ~ dmnorm(0,0.0001)
log.mu.k ~ dunif(0,1)
mu.t0 ~ dmnorm(0,0.0001)

# - Put on natural scale for reporting
mu.Linf = exp(log.mu.Linf)
mu.k = exp(log.mu.k)
 
}"

# write model to a text file
writeLines(jags.mod, "jags_model.txt")

##### PARAMETERS TO MONITOR #####
params = c("mu.Linf", "mu.k", "mu.t0", "Linf.group", "k.group", "t0.group","Linf.i", "k.i", "t0.i", 
           "ind.cov.mat","group.cov.mat","sigma", 
           "log.linf.beta.river", "log.k.beta.river", "t0.beta.river", "log.linf.beta.sex", "log.k.beta.sex", "t0.beta.sex" )


#### A#ssign data to list ##### 
dat = list(age = age, length = length, sex = sex, river = river, sample.id = sample.id, N = N, G = G, Nind = Nind, group = group, z.group.mat = diag(x=1,nrow=3), z.ind.mat = diag(x=1,nrow=3), mu.re = rep(0, 3))


##### MCMC DIMENSIONS #####
ni = 50000
nb = 100000
na = 100000
nt = 1000
nc = 4
n.iter = ni + nb

##### RUN THE MODEL IN JAGS #####
runJagsOut <- run.jags( model="jags_model.txt" ,
                        monitor=params ,
                        data=dat ,
                        n.chains=nc ,
                        adapt=na ,
                        burnin=nb ,
                        sample=ni ,
                        thin=nt ,
                        method = "parallel",
                        plots=FALSE )

plot(runJagsOut) # Converged?

summ <- summary(runJagsOut)