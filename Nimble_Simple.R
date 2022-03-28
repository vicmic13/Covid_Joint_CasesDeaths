# Write the code for the joint model p(cases, deaths)=p(cases)*p(deaths | cases)
Simple_Code <- nimbleCode({
  # Likelihood
  for(i in 1:N){
    cases[i] ~ dpois(lambda_cases[i])
    log(lambda_cases[i]) <- log(E_cases[i]) + beta0_cases + inprod(X[i,1:p],beta_cases[1:p])
    censored[i] ~ dinterval(deaths[i], c[i,])
    deaths[i] ~ dpois(lambda_deaths[i])
    log(lambda_deaths[i]) <- log(E_deaths[i]) + beta0_deaths + inprod(X[i,1:p],beta_deaths[1:p])
  }
  
  # Priors
  beta0_cases ~ dnorm(0, sd=10)
  beta0_deaths ~ dnorm(0, sd=10)
  for(k in 1:p){
    beta_cases[k] ~ dnorm(0, sd=10)
    beta_deaths[k] ~ dnorm(0, sd=10)
  }
})

# Give Nimble all the data information
Consts_Simple = list(N=N,
                     E_cases=data$E_cases,
                     E_deaths=data$E_cases_deaths,
                     X=matrix(c(data$beds_log, data$diploma, data$age), N, 3),
                     p=3,
                     c=c)
Data_Simple = list(cases = data$CASES, deaths=data$deaths, censored=data$cens)
Inits_Simple = list(beta0_deaths=0, beta_deaths=rep(0,3), deaths=rep(1,N),
                    beta0_cases=0, beta_cases=rep(0,3))
