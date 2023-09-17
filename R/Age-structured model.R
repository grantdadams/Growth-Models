## Load packages ----
# install.packages("pacman")
pacman::p_load(forecast)


## Define population and observation model ----
# Mortality and population specifications
M <- 0.2
nyrs <- 50
yrs <- 1:nyrs
nages <- 10
ages <- 1:nages

# - Weight at age
winf <- 10
vb_k <- 0.25
wt <- winf * (1-exp(-vb_k * ages))

# - Maturity
mat <- as.numeric(ages>4) # Knife-edge maturity above age 4

# - Recruitment
Rhat <- 1000 
Rdev <- arima.sim(list(order=c(1,0,0), ar=.5), n = nyrs, rand.gen = rlnorm) # Follows AR1 process
R <- Rhat * Rdev
plot(y = R, x = yrs, type = "l")

# - Fishing mortality and selectivity
# -- here I am assuming the seperability where there is an increasing and decreasing fishing mortality
# -- and logistic fishery selectivity
fsh_mort <- c(seq(0, 0.2, length.out = nyrs-10), seq(0.2, 0.1, length.out = 10)) # Increasing and decreasing F
fsh_slope <- 2
fsh_asymp <- 4
fsh_sel <- 1/(1+exp(-fsh_slope * (ages-fsh_asymp)))
plot(x = ages, y = fsh_sel, type = "l")

# - Survey selectivity and catchability
srv_q <- 0.4
srv_slope <- 2
srv_asymp <- 3
srv_sel <- 1/(1+exp(-srv_slope * (ages-srv_asymp))) # Logistic selectivity that targets young fish
plot(x = ages, y = srv_sel, type = "l")


## Fill population and observation model ----
n_at_age <- matrix(0, nrow = nages, ncol = nyrs + 1)
colnames(n_at_age) <- 0:nyrs; rownames(n_at_age) <- ages
srv_comp <- fsh_comp <- matrix(0, nrow = nages, ncol = nyrs)
biomass <- srv_biom <- ssb <- catch <- c()

# - Calculate initial abundance and recruitment
n_at_age[1,2:(nyrs+1)] <- R # Age-1 is recruits
n_at_age[1,1] <- Rhat # Age-1, Yr-0 is Rhat
n_at_age[2:nages,1] <- Rhat * exp(-M * 1:(nages-1))
n_at_age[nages,1] <- n_at_age[nages,1]/(1-exp(-M*nages)) # Max-age, solve for geometric series for equilibrium

# - Loop through years and ages
for(y in 2:nyrs){ # Year one starts at equilibrium
  for(a in 2:nages){
    n_at_age[a,y] <- n_at_age[a-1,y-1] * exp(-M - fsh_sel[a-1] * fsh_mort[y-1])
  }
  
  # - Calculate biomass
  biomass[y] <- sum(n_at_age[,y] * wt)
  ssb[y] <- sum(n_at_age[,y] * wt * mat * 0.5) # Assuming 50% females 
  
  
  # - Calculate fishery and survey age composition
  srv_comp[,y] <- n_at_age[,y] * srv_sel # Assuming happens at month-0
  fsh_comp[,y] <- n_at_age[,y] * (fsh_sel * Fmort[y])/ (M + fsh_sel * fsh_mort[y]) * (1-exp(-M - fsh_sel * fsh_mort[y])) # Baranov catch equation
  # - Baranov catch equation is the integral of dN/dT -(F + M) times the ratio of F to (M+F)
  # - The 1 - exp(-M+F) is calculating the number of fish lost due to F + M: (N-N*exp(-M+F))
  
  # - Calculate catch and survey biomass
  srv_biom[y] <- sum(srv_comp[,y] * wt * srv_q) 
  catch[y] <- sum(fsh_comp[,y] * wt) 
  
  # - Normalize comp data
  srv_comp[,y] <- srv_comp[,y]/sum(srv_comp[,y])
  fsh_comp[,y] <- fsh_comp[,y]/sum(fsh_comp[,y])
}


## Add observation error ----
# - Simulate observed catch and survey data from lognormal distribution
srv_obs <- rlnorm(n = nyrs, meanlog = log(srv_biom), sd = 0.2)
catch_obs <- rlnorm(n = nyrs, meanlog = log(catch), sd = 0.1)

# - Simulate age composition data from multinomial with assumed sample size of 200 for survey and 100 for fishery
# -- Size is the multinomial sample size (See Hamel and Stewart 2014 for calculation with observed data)
srv_comp_obs <- apply(srv_comp, 1, function(x) rmultinom(n = 1, size = 200, prob = x)) 
srv_comp_obs <- apply(srv_comp_obs, 1, function(x) x/sum(x)) # Normalize

fsh_comp_obs <- apply(fsh_comp, 1, function(x) rmultinom(n = 1, size = 100, prob = x))
fsh_comp_obs <- apply(fsh_comp_obs, 1, function(x) x/sum(x))
