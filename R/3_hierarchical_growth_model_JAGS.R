####################################################################################
# Simulate length-at-age data ######################################################
####################################################################################
# Set Seed
set.seed(1234)

# Number of Groups
Groups=6

# True VBGM hyperparameters c(mean, sd)
true.Linf = c(500,20)
true.k = c(.3,.05)
true.t0 = c(1.5,.4)

sigma <- 15

# Age range
ages = seq(from=1,to=15, by = .05)

# Empty matrix and vectors to fill with parameters and data, respectively
param.mat = matrix(NA,Groups,3,byrow = T)
colnames(param.mat) <- c("Linf", "k", "t0")
ctr = 0
age = c()
length = c()
group = c()

# Simulate group level parameters
for(i in 1:Groups){
  param.mat[i,1] = rnorm(1,true.Linf[1],true.Linf[2]) # Assign group level Linf
  param.mat[i,2] = rnorm(1,true.k[1],true.k[2]) # Assign group level k
  param.mat[i,3] = rnorm(1,true.t0[1],true.t0[2]) # Assign group level t0
  n.samples = sample(200:1000, 1) # Number of samples per group s
  
  # Simulate data
  for(j in 1:n.samples) {
    ctr = ctr + 1 # Indexing variable
    age[ctr] = sample(ages, 1) # Sample randon age from age range
    length[ctr] = (param.mat[i,1] * (1 - exp(-param.mat[i,2]*(age[ctr]-param.mat[i,3])))) + rnorm(1,0,sigma)
    group[ctr] = i
  }
}

# Assign data to data frame
dat = data.frame(age = age, length = length, group = group, N = length(age), G = length(unique(group)))
dat <- dat[which(dat$length > 0),]

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
tau.y = pow(sigma, -2)
sigma ~ dunif(0,100)
 
# Level-2 parameters
for(j in 1:G){
Linf[j] ~ dnorm(mu.Linf, tau.Linf)
k[j] ~ dnorm(mu.k, tau.k)
t0[j] ~ dnorm(mu.t0, tau.t0)
}
 
# Priors for level-2 parameters
log.mu.Linf ~ dnorm(0,0.0001)
log.mu.k ~ dnorm(0,0.0001)
mu.t0 ~ dnorm(0,0.0001)
 
# Get hyperparameters on untransformed scale
mu.Linf = exp(log.mu.Linf)
mu.k = exp(log.mu.k)
 
# Precision
tau.Linf = pow(sig.Linf,-2)
tau.k = pow(sig.k,-2)
tau.t0 = pow(sig.t0,-2)
 
# SD of parameters
sig.Linf ~ dunif(0,10)
sig.k ~ dunif(0,10)
sig.t0 ~ dunif(0,10)
}"

# write model to a text file
writeLines(jags.mod, "jags_model.txt")

##### PARAMETERS TO MONITOR #####
params = c("Linf", "k", "t0", "mu.Linf", "mu.k", "mu.t0", "mu.Linf", "mu.k", "mu.t0", "sig.Linf","sig.k","sig.t0","sigma" )

##### MCMC DIMENSIONS #####
ni = 5000
nb = 5000
na = 5000
nt = 100
nc = 3
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
                        summarise=FALSE ,
                        plots=FALSE )

plot(runJagsOut) # Converged?

summ <- summary(runJagsOut)
summ[,1:4]
param.mat
true.Linf
true.k
true.t0
