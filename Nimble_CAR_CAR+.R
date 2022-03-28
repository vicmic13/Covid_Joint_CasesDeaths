# Write the code for the joint model p(cases, deaths)=p(cases)*p(deaths | cases)
CAR_Code <- nimbleCode({
  for(i in 1:N){
    # Likelihood
    # Cases
    cases[i] ~ dpois(exp(log(E_cases[i]) + beta0_cases + inprod(X[i,1:p], beta_cases[1:p]) + s[i,1]))
    lambda_cases[i] <- exp(log(E_cases[i]) + beta0_cases + inprod(X[i,1:p], beta_cases[1:p]) + s[i,1])
    # CAR model: X is the matrix of covariates in their original scale
    # CAR+ model: X is the matrix of covariates after adjusting for spatial confounding (Spatial+)
    # i.e., for CAR+, X_ik = r_ik, following the paper notation
    
    # Deaths
    censored[i] ~ dinterval(deaths[i], c[i,])
    deaths[i] ~ dpois(exp(log(E_deaths[i]) + beta0_deaths + inprod(X[i,1:p], beta_deaths[1:p]) + s[i,2]))
    lambda_deaths[i] <- exp(log(E_deaths[i]) + beta0_deaths + inprod(X[i,1:p], beta_deaths[1:p]) + s[i,2])
    # CAR model: X is the matrix of covariates in their original scale
    # CAR+ model: X is the matrix of covariates after adjusting for spatial confounding (Spatial+)
    # i.e., for CAR+, X_ik = r_ik, following the paper notation
    
    # Joint latent effects
    s[i,1:2] <- Achol[1:2,1:2] %*% u[i, 1:2]
    s_cases[i] <- s[i,1]
    s_deaths[i] <- s[i,2]
  }
  
  # Priors
  beta0_cases ~ dnorm(0, sd=10)
  beta0_deaths ~ dnorm(0, sd=10)
  
  # Spatial effects independently defined for cases and deaths
  u[1:N,1] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], 1, zero_mean=1)
  u[1:N,2] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], 1, zero_mean=1)
  # Spatial structure defined through vectors adj, weights and num
  # L = number of pairs of neighbours; weights = vector of 1s
  # adj = all the neighbours for each area; num = number of neighbours for each area
  # e.g. num[1]=2: area 1 has 2 neighbours, which are adj[1] and adj[2]
  #      num[2]=2: area 2 has 2 neighbours, which are adj[3] and adj[4], ...
  
  # Define latent effects' covariance matrix
  Sigma_b[1,1] <- sigma_b_cases^2
  Sigma_b[1,2] <- rho*sigma_b_cases*sigma_b_deaths
  Sigma_b[2,1] <- rho*sigma_b_cases*sigma_b_deaths
  Sigma_b[2,2] <- sigma_b_deaths^2
  Achol[1:2,1:2] <- t(chol(Sigma_b[1:2,1:2]))
  
  sigma_b_cases ~ T(dt(0,1,1),0,)
  sigma_b_deaths ~ T(dt(0,1,1),0,)
  rho ~ dunif(-1,1)
  
  for(k in 1:p){
    beta_cases[k] ~ dnorm(0, sd=10)
    beta_deaths[k] ~ dnorm(0, sd=10)
  }
})

# Give Nimble all the data information
# Constants for the CAR model
Consts_CAR = list(N=N,
                  E_cases=data$E_cases,
                  E_deaths=data$E_cases_deaths,
                  X=matrix(c(data$beds_log, data$diploma, data$age), N, 3),
                  p=3,
                  L=length(neigh),
                  adj=neigh,
                  weights=rep(1,length(neigh)),
                  num=numneigh,
                  c=c)

# Constants for the CAR+ model
Consts_CAR_r = list(N=N,
                    E_cases=data$E_cases,
                    E_deaths=data$E_cases_deaths,
                    X=r_x_pm,
                    p=3,
                    L=length(neigh),
                    adj=neigh,
                    weights=rep(1,length(neigh)),
                    num=numneigh,
                    c=c)

Data_CAR = list(cases = data$CASES, deaths=data$deaths, censored=data$cens)
Inits_CAR = list(beta0_deaths=0, beta_deaths=rep(0,3), deaths=rep(1,N),
                 beta0_cases=0, beta_cases=rep(0,3), 
                 s_deaths=rep(0,N), s_cases=rep(0,N),
                 u=matrix(0,nrow=N,ncol=2),
                 sigma_b_cases=1, sigma_b_deaths=1, rho=0)
