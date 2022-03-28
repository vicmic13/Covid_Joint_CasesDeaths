setwd(".")
##########################
#  PACKAGES & FUNCTIONS  #
##########################
library(ggplot2)
library(ggpubr)
library(rgdal)
library(broom) # for tidy: extract geo data
library(dplyr) # for left_join: add data to geo data
library(spdep)
library(nimble)
library(coda)

adjlist = function(W,N){ 
  adj=0
  for(i in 1:N){  
    for(j in 1:N){  
      if(W[i,j]==1){adj = append(adj,j)}
    }
  }
  adj = adj[-1]
  return(adj)}

###########################
#         DATASETS        #
###########################
##### Shapefiles
Mtl_sh = readOGR("./Datasets", "LIMADMIN", stringsAsFactors = FALSE)
Mtl_neigh=tidy(Mtl_sh)  
Mtl_sh$id=row.names(Mtl_sh) ; Mtl_neigh=left_join(Mtl_neigh, Mtl_sh@data, by="id")

# Remove Ile-Dorval
Mtl_neigh = Mtl_neigh[Mtl_neigh$NOM != "L'Ile-Dorval",]
N = length(unique(Mtl_neigh$NOM))

##### Dataset
data = read.csv("./Datasets/Data.csv", header=T, sep=",")
data$E_cases = data$POPULATION*(sum(data$CASES)/sum(data$POPULATION))
data$income = as.numeric(scale(data$INCOME))
data$diploma = (data$DIPLOMA)/100
data$age = as.numeric(scale(data$AGE))
data$deaths = as.numeric(data$DEATHS)
data$deaths_noNA=data$deaths ; data[is.na(data$deaths),"deaths_noNA"] = c(0,0)
data$E_cases_deaths = data$CASES*(sum(data$deaths_noNA)/sum(data$CASES))
data$beds = as.numeric(scale(data$BEDS))
data$beds_log=as.numeric(scale(log(data$BEDS+1)))
data$cens = rep(1, N)
c=cbind(c(rep(-1,2), 0, rep(-1,16), 0, rep(-1,13)), c(rep(Inf,2), 4, rep(Inf,16), 4, rep(Inf,13)))

##### Adjacency matrix
W = read.csv("./Datasets/Adjacency.csv", header=T, row.names=1, sep=",")
colnames(W) <- rownames(W)

# Spatial structure to use in Nimble
neigh = adjlist(W, N)
numneigh = apply(W,2,sum)

# Merge the datasets
Mtl_neigh = left_join(Mtl_neigh[,c(1,2,6,11)], data[,-c(3,5,6,7,9)], by="NOM")

# Some map
ggplot() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "bottom") + 
  geom_polygon(data=Mtl_neigh, aes(x = long, y = lat, group = group, fill=age), col="black")+
  scale_fill_gradient("Age", low="white", high="gray37") +
  coord_equal()

#######################
#    NIMBLE: Simple   #
#######################
source("./Nimble_Simple.R")
Simple_Model = nimbleModel(code=Simple_Code, 
                           constants=Consts_Simple, 
                           data=Data_Simple, 
                           inits=Inits_Simple)

# Define which parameters to monitor (e.g. to extract posterior samples for)
Simple_conf <- configureMCMC(Simple_Model, monitors = c("beta0_deaths", "beta_deaths","deaths",
                                                        "beta0_cases", "beta_cases", "cases",
                                                        "lambda_cases", "lambda_deaths"))
Simple_MCMC=buildMCMC(Simple_conf, enableWAIC=T)
Simple_Cmodel=compileNimble(Simple_Model)
Simple_Cmcmc=compileNimble(Simple_MCMC, project = Simple_Model)
Simple_runMCMC=runMCMC(Simple_Cmcmc, nburnin=25000, niter=50000, thin=10, nchains=2, samplesAsCodaMCMC = TRUE, WAIC=T)

# Check convergence:
plot(Simple_runMCMC$samples[,c("lambda_deaths[20]")])

#######################
#     NIMBLE: IID     #
#######################
source("./Nimble_IID.R")

IID_Model = nimbleModel(code=IID_Code, 
                        constants=Consts_IID, 
                        data=Data_IID, 
                        inits=Inits_IID)

# Define which parameters to monitor (e.g. to extract posterior samples for)
IID_conf <- configureMCMC(IID_Model, monitors = c("beta0_deaths", "beta_deaths",
                                                  "b_deaths", "sigma_b_deaths",  "deaths",
                                                  "beta0_cases", "beta_cases",
                                                  "lambda_cases", "lambda_deaths",
                                                  "b_cases", "sigma_b_cases", "cases", "rho"))

IID_MCMC=buildMCMC(IID_conf, enableWAIC=T)
IID_Cmodel=compileNimble(IID_Model)
IID_Cmcmc=compileNimble(IID_MCMC, project = IID_Model)
IID_runMCMC=runMCMC(IID_Cmcmc, nburnin=5000, niter=20000, thin=5, nchains=2, samplesAsCodaMCMC = TRUE, WAIC=T)

# Check convergence:
plot(IID_runMCMC$samples[,c("beta0_deaths")])

#######################
#   NIMBLE: Spatial+  #
#######################
source("./Nimble_Spatial+.R")

SpatialX_Model = nimbleModel(code=SpatialX_Code, 
                             constants=Consts_SpatialX, 
                             data=Data_SpatialX, 
                             inits=Inits_SpatialX)

# Define which parameters to monitor (e.g. to extract posterior samples for)
SpatialX_conf <- configureMCMC(SpatialX_Model, monitors = c("beta0_x", "f_x",
                                                            "sigma_e", "sigma_f", "r"))

SpatialX_MCMC=buildMCMC(SpatialX_conf, enableWAIC=T)
SpatialX_Cmodel=compileNimble(SpatialX_Model)
SpatialX_Cmcmc=compileNimble(SpatialX_MCMC, project = SpatialX_Model)
SpatialX_runMCMC=runMCMC(SpatialX_Cmcmc, nburnin=10000, niter=20000, thin=10, nchains=2, samplesAsCodaMCMC = TRUE, WAIC=T)

# Check convergence:
plot(SpatialX_runMCMC$samples[,c("r[1, 1]")])

# Get residuals' posterior means to replace as covariates in CAR+ and BYM+ models
Samples_SpatialX = rbind(SpatialX_runMCMC$samples[[1]], SpatialX_runMCMC$samples[[2]])
r_x1 = Samples_SpatialX[,103:135]
r_x2 = Samples_SpatialX[,136:168]
r_x3 = Samples_SpatialX[,169:201]

r_x_pm = cbind(apply(r_x1,2,mean), apply(r_x2,2,mean), apply(r_x3,2,mean))
colnames(r_x_pm)=c("res_beds_log","res_diploma","res_age")

#######################
#  NIMBLE: CAR & CAR+ #
#######################
source("./Nimble_CAR_CAR+.R")

##### RUN CAR
CAR_Model = nimbleModel(code=CAR_Code, 
                        constants=Consts_CAR, 
                        data=Data_CAR, 
                        inits=Inits_CAR)
CAR_conf <- configureMCMC(CAR_Model, monitors = c("beta0_deaths", "beta_deaths", "s_deaths", "deaths",
                                                  "beta0_cases", "beta_cases", "s_cases","u",
                                                  "sigma_b_cases", "sigma_b_deaths", "rho",
                                                  "cases", "lambda_cases", "lambda_deaths"))

# Change some of the default samplers to gain convergence efficiency 
CAR_conf$removeSampler(c('beta0_cases', "beta0_deaths", "beta_cases", "beta_deaths"))
CAR_conf$addSampler(target = c('beta0_cases'), type = "ess")
CAR_conf$addSampler(target = c('beta0_deaths'), type = "ess")
CAR_conf$addSampler(target = c('beta_cases[1]'), type = "ess")
CAR_conf$addSampler(target = c('beta_cases[2]'), type = "ess")
CAR_conf$addSampler(target = c('beta_cases[3]'), type = "ess")
CAR_conf$addSampler(target = c('beta_deaths[1]'), type = "ess")
CAR_conf$addSampler(target = c('beta_deaths[2]'), type = "ess")
CAR_conf$addSampler(target = c('beta_deaths[3]'), type = "ess")
CAR_MCMC=buildMCMC(CAR_conf, enableWAIC=T)
CAR_Cmodel=compileNimble(CAR_Model)
CAR_Cmcmc=compileNimble(CAR_MCMC, project = CAR_Model)
CAR_runMCMC=runMCMC(CAR_Cmcmc, nburnin=5000, niter=20000, thin=5, nchains=2, samplesAsCodaMCMC = TRUE, WAIC=T)

# Check convergence:
plot(CAR_runMCMC$samples[,c("beta_deaths[2]")])

##### RUN CAR+
CAR_Model_r = nimbleModel(code=CAR_Code, 
                          constants=Consts_CAR_r, 
                          data=Data_CAR, 
                          inits=Inits_CAR)
CAR_conf_r <- configureMCMC(CAR_Model_r, monitors = c("beta0_deaths", "beta_deaths", "s_deaths", "deaths",
                                                      "beta0_cases", "beta_cases", "s_cases", "u",
                                                      "sigma_b_cases", "sigma_b_deaths", "rho",
                                                      "cases", "lambda_cases", "lambda_deaths"))
CAR_MCMC_r=buildMCMC(CAR_conf_r, enableWAIC=T)
CAR_Cmodel_r=compileNimble(CAR_Model_r)
CAR_Cmcmc_r=compileNimble(CAR_MCMC_r, project = CAR_Model_r)
CAR_runMCMC_r=runMCMC(CAR_Cmcmc_r, nburnin=400000, niter=800000, thin=160, nchains=2, samplesAsCodaMCMC = TRUE, WAIC=T)

# Check convergence:
plot(CAR_runMCMC_r$samples[,c("beta_deaths[2]")])

#######################
#  NIMBLE: BYM & BYM+ #
#######################
source("./Nimble_BYM_BYM+.R")

##### RUN BYM
BYM_Model = nimbleModel(code=BYM_Code, 
                        constants=Consts_BYM, 
                        data=Data_BYM, 
                        inits=Inits_BYM)
BYM_conf <- configureMCMC(BYM_Model, monitors = c("beta0_deaths", "beta_deaths", "s_deaths", "theta_deaths",  
                                                  "sigma_theta_deaths", "sigma_s_deaths", "deaths",
                                                  "beta0_cases", "beta_cases", "s_cases", "theta_cases", 
                                                  "sigma_theta_cases", "sigma_s_cases", "cases",
                                                  "lambda_cases", "lambda_deaths",
                                                  "rho_s", "rho_theta"))
BYM_MCMC=buildMCMC(BYM_conf, enableWAIC=T)
BYM_Cmodel=compileNimble(BYM_Model)
BYM_Cmcmc=compileNimble(BYM_MCMC, project = BYM_Model)
BYM_runMCMC=runMCMC(BYM_Cmcmc, nburnin=5000, niter=20000, thin=5, nchains=2, samplesAsCodaMCMC = TRUE, WAIC=T)

# Check convergence:
plot(BYM_runMCMC$samples[,c("beta_deaths[2]")])

##### RUN BYM+
BYM_Model_r = nimbleModel(code=BYM_Code, 
                          constants=Consts_BYM_r, 
                          data=Data_BYM, 
                          inits=Inits_BYM)
BYM_conf_r <- configureMCMC(BYM_Model_r, monitors = c("beta0_deaths", "beta_deaths", "s_deaths", "theta_deaths",  
                                                      "sigma_theta_deaths", "sigma_s_deaths", "deaths",
                                                      "beta0_cases", "beta_cases", "s_cases", "theta_cases", 
                                                      "sigma_theta_cases", "sigma_s_cases", "cases",
                                                      "lambda_cases", "lambda_deaths",
                                                      "rho_s", "rho_theta"))
BYM_MCMC_r=buildMCMC(BYM_conf_r, enableWAIC=T)
BYM_Cmodel_r=compileNimble(BYM_Model_r)
BYM_Cmcmc_r=compileNimble(BYM_MCMC_r, project = BYM_Model_r)
BYM_runMCMC_r=runMCMC(BYM_Cmcmc_r, nburnin=400000, niter=800000, thin=160, nchains=2, samplesAsCodaMCMC = TRUE, WAIC=T)

# Check convergence:
plot(BYM_runMCMC_r$samples[,c("beta_deaths[2]")])

#######################
#       RESULTS       #
#######################
samples_Simple = rbind(Simple_runMCMC$samples[[1]], Simple_runMCMC$samples[[2]])
samples_IID = rbind(IID_runMCMC$samples[[1]], IID_runMCMC$samples[[2]])
samples_CAR = rbind(CAR_runMCMC$samples[[1]], CAR_runMCMC$samples[[2]])
samples_BYM = rbind(BYM_runMCMC$samples[[1]], BYM_runMCMC$samples[[2]])
samples_CAR_r = rbind(CAR_runMCMC_r$samples[[1]], CAR_runMCMC_r$samples[[2]])
samples_BYM_r = rbind(BYM_runMCMC_r$samples[[1]], BYM_runMCMC_r$samples[[2]])


