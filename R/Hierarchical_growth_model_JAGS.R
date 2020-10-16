####################################################################################
# Simulate length-at-age data ######################################################
####################################################################################
# Set Seed
set.seed(1234)

# Number of Groups
Groups=6

# True VBGM parameters c(mean, sd)
true.Linf = c(500,20)
true.k = c(.3,.05)
true.t0 = c(1.5,.4)

sigma

# Age range
ages = seq(from=1,to=15, by = .05)

# Empty matrix and vectors to fill with parameters and data, respectively
param.mat = matrix(NA,6,3,byrow = T)
ctr = 0
age = c()
length = c()
group = c()

# Simulate Data
for(i in 1:Groups){
  param.mat[i,1] = rnorm(1,true.Linf[1],true.Linf[2]) # Assign group level Linf
  param.mat[i,2] = rnorm(1,true.k[1],true.k[2]) # Assign group level k
  param.mat[i,3] = rnorm(1,true.t0[1],true.t0[2]) # Assign group level t0
  n.samples = sample(200:1000, 1) # Number of samples per group s
  for(j in 1:n.samples) {
    ctr = ctr + 1 # Indexing variable
    age[ctr] = sample(ages, 1) # Sample randon age from age range
    length[ctr] = (param.mat[i,1] * (1 - exp(-param.mat[i,2]*(age[ctr]-param.mat[i,3])))) + rnorm(1,0,sigma)
    group[ctr] = i
  }
}

# Assign data to data frame
dat = data.frame(age = age, length = length, group = group, N = length(age), G = length(unique(group)))

# Plot the data
for (i in 1:length(unique(dat$group))){
  sub = dat[which(dat$group==i),]
  if ( i == 1){
    plot(sub$age, sub$length, col = i, xlab = "Age", ylab = "Length")
  } else {
    points(sub$age, sub$length , col = i)
  }
}

# Assign data to list
dat = list(age = age, length = length, group = group, N = length(age), G = length(unique(group)))

#################################################################################### 
# JAGS VBGM ######################################################################## 
####################################################################################
require(runjags)
require(coda)

##### SPECIFY JAGS MODEL CODE #####
jags.mod =
  "model{
for(i in 1:N){
length[i] ~ dnorm(y.hat[i], tau.y)
y.hat[i] = Linf[group[i]] * (1-exp(-k[group[i]] * (age[i] - t0[group[i]] )))
}
 
# SD
tau.y sigma ~ dunif(0,100)
 
# Level-2 parameters
for(j in 1:G){
Linf[j] k[j] t0[j] B.raw.hat[j,1] ~ dnorm(mu.Linf.raw, tau.Linf)
B.raw.hat[j,2] ~ dnorm(mu.k.raw, tau.k)
B.raw.hat[j,3] ~ dnorm(mu.t0.raw, tau.t0)
}
 
#priors for level-2 parameters
mu.Linf.raw ~ dnorm(0,0.0001)
mu.k.raw ~ dnorm(0,0.0001)
mu.t0.raw ~ dnorm(0,0.0001)
 
# Get hyperparameters on untransformed scale
mu.Linf mu.k mu.t0
 
# Precision
tau.Linf = pow(sig.Linf,-2)
tau.k = pow(sig.Linf,-2)
tau.t0 = pow(sig.Linf,-2)
 
# SD of parameters
sig.Linf ~ dunif(0,10)
sig.k ~ dunif(0,10)
sig.t0 ~ dunif(0,10)
}"

# write model to a text file
writeLines(jags.mod, "jags_model.txt")

##### PARAMETERS TO MONITOR #####
params = c("Linf", "k", "t0", "mu.Linf", "mu.k", "mu.t0", "mu.Linf.raw", "mu.k.raw", "mu.t0.raw", "sig.Linf","sig.k","sig.t0","sigma" )

##### MCMC DIMENSIONS #####
ni = 500
nb = 2000
na = 1000
nt = 10
nc = 3
n.iter = ni + nb

##### RUN THE MODEL IN JAGS #####
runJagsOut model="jags_model.txt" ,
monitor=params ,
data=dat ,
n.chains=nc ,
adapt=na ,
burnin=nb ,
sample=ni ,
thin=nt ,
summarise=FALSE ,
plots=FALSE )
summary(runJagsOut)