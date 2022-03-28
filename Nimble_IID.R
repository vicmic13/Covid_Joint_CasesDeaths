# Write the code for the joint model p(cases, deaths)=p(cases)*p(deaths | cases)
IID_Code <- nimbleCode({
  for(i in 1:N){
    # Likelihood
    # Cases
    cases[i] ~ dpois(exp(b[i,1]))
    lambda_cases[i] <- exp(b[i,1])
    
    # Deaths
    censored[i] ~ dinterval(deaths[i], c[i,])
    deaths[i] ~ dpois(exp(b[i,2]))
    lambda_deaths[i] <- exp(b[i,2])
    
    # Priors
    # Joint latent effects
    b[i,1:2] ~ dmnorm(mean=mean_b[i,1:2], cholesky=Chol_Sigma_b[1:2,1:2], prec_param=0)
    b_cases[i] <- b[i,1] - mean_b[i,1]
    b_deaths[i] <- b[i,2] - mean_b[i,2]
    
    # Nimble converges faster with centered structures
    mean_b[i,1] <- log(E_cases[i]) + beta0_cases + inprod(X[i,1:p], beta_cases[1:p])
    mean_b[i,2] <- log(E_deaths[i]) + beta0_deaths + inprod(X[i,1:p], beta_deaths[1:p])
  }
  
  # Define latent effects' covariance matrix
  Sigma_b[1,1] <- sigma_b_cases^2
  Sigma_b[1,2] <- rho*sigma_b_cases*sigma_b_deaths
  Sigma_b[2,1] <- rho*sigma_b_cases*sigma_b_deaths
  Sigma_b[2,2] <- sigma_b_deaths^2
  Chol_Sigma_b[1:2,1:2] <- chol(Sigma_b[1:2,1:2])
  
  # All the other parameters
  beta0_cases ~ dnorm(0, sd=10)
  beta0_deaths ~ dnorm(0, sd=10)
  for(k in 1:p){
    beta_cases[k] ~ dnorm(0, sd=10)
    beta_deaths[k] ~ dnorm(0, sd=10)
  }
  sigma_b_cases ~ T(dt(0,1,1),0,)
  sigma_b_deaths ~ T(dt(0,1,1),0,)
  rho ~ dunif(-1,1)
})

# Give Nimble all the data information
Consts_IID = list(N=N,
                  E_cases=data$E_cases,
                  E_deaths=data$E_cases_deaths,
                  X=matrix(c(data$beds_log, data$diploma, data$age), N, 3),
                  p=3,
                  c=c)
Data_IID = list(cases = data$CASES, deaths=data$deaths, censored=data$cens)
Inits_IID = list(beta0_deaths=0, beta_deaths=rep(0,3), deaths=rep(1,N),
                 beta0_cases=0, beta_cases=rep(0,3), 
                 b=matrix(0, nrow=N, ncol=2), sigma_b_cases=1, sigma_b_deaths=1, rho=0)
