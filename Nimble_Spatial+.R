# Write the code for the Spatial+ model
SpatialX_Code <- nimbleCode({
  for(j in 1:K){
    # K is the number of covariates, 3 here
    for(i in 1:N){
      X[i,j] ~ dnorm(beta0_x[j] + f_x[i,j], sd = sigma_e[j])
      
      # We want to take the residuals in our CAR+ and BYM+ models, so monitor the r's:
      r[i,j] <- X[i,j] - beta0_x[j] - f_x[i,j]
    }
    f_x[1:N,j] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau[j], zero_mean=1)   
    tau[j] <- 1/(sigma_f[j])^2
    # Spatial structure defined through vectors adj, weights and num
    # L = number of pairs of neighbours; weights = vector of 1's
    # adj = all the neighbours for each area; num = number of neighbours for each area
    # e.g. num[1]=2: area 1 has 2 neighbours, which are adj[1] and adj[2]
    # num[2]=2: area 2 has 2 neighbours, which are adj[4] and adj[5]...
    
    # Priors
    beta0_x[j] ~ dnorm(0, sd=10)
    sigma_f[j] ~ T(dt(0,1,1),0,)
    sigma_e[j] ~ T(dt(0,1,1),0,)
  }
})

# Give Nimble all the data information
Consts_SpatialX = list(K=3,
                       N=N,
                       L=length(neigh),
                       adj=neigh,
                       weights=rep(1,length(neigh)),
                       num=numneigh)
Data_SpatialX = list(X=matrix(c(data$beds_log, data$diploma, data$age), N, 3))
Inits_SpatialX = list(beta0_x=rep(0,3), f_x=matrix(0, N, 3), sigma_e=rep(1,3), sigma_f=rep(1,3))
