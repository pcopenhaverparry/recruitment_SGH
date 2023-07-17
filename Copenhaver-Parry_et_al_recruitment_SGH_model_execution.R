#Model Execution Code
#Copenhaver-Parry, P.E., Byerly, S.N., Greenwood, J.W., Retter, C.J., Beedlow, 
#P.A., & Lee, E.H. Simultaneous changes in climate and competition limit forest
#recruitment with increasing climatic stress.

rm(list=ls())

#Load required packages
library(rstan)
library(loo)

#Select species to model
#species <- 'TSHE'
species <- 'PSME'

#Specify model file
model_file <- 'Copenhaver-Parry_et_al_recruitment_SGH_model.stan'

####Prepare data####

#Load data
data <- read.csv(paste('Copenhaver-Parry_et_al_recruitment_SGH_',species, '_data.csv',sep=''))

data <- na.omit(data)

#Extract covariates from data and assign to covariate matrix 
X.unscaled <- cbind(C=data[,'Veg_Cover'],
                    GSF=data[,'GSF'],
                    BA=data[,'BA'],
                    Germ_AT=data[,'Germ_AT'],
                    Germ_PPT=data[,'Germ_PPT'],
                    Germ_ATmax=data[,'Germ_ATmax'],
                    Germ_ATmin=data[,'Germ_ATmin'],
                    Winter_AT=data[,'Winter_AT'],
                    Winter_PPT=data[,'Winter_PPT'],
                    AprJun_AT=data[,'AprJun_AT'],
                    AprJun_PPT=data[,'AprJun_PPT'],
                    Summer_ATmax=data[,'Summer_ATmax'],
                    Summer_PPT=data[,'Summer_PPT'])

#Standardize covariates
means <- vector()
sds <- vector()
for(i in 1:ncol(X.unscaled)){
  means[i] <- mean(X.unscaled[,i])
  sds[i] <- sd(X.unscaled[,i])
}
means.sds <- rbind(means,sds)
colnames(means.sds) <- colnames(X.unscaled)

X <- matrix(data=NA, nrow=nrow(X.unscaled), ncol=ncol(X.unscaled))
colnames(X) <- colnames(X.unscaled)

for(i in 1:ncol(X.unscaled)){
  X[,i] <- (X.unscaled[,i]-means.sds[1,i])/means.sds[2,i]
}

#Interaction terms
intX <- expand.grid(colnames(X),colnames(X))

for(j in 1:nrow(intX)){
  xnames <- colnames(X)
  X <- cbind(X, X[,intX[j,1]]*X[,intX[j,2]])
  colnames(X)<- c(xnames,paste(intX[j,1],intX[j,2],sep="x"))
}

X <- cbind(X,X2)

####Set up candidate models####

#Load table of candidate models
mods.to.test <- read.csv('Copenhaver-Parry_et_al_recruitment_SGH_mods.csv')

#Organize candidate models
X.keep<-colnames(mods.to.test)

names.keep <- colnames(X)[colnames(X) %in% X.keep] 

x<-X[,names.keep]

model_names <- c('germination', 'drought','frost', 'chilling','seeds','summer')

####Fit models####

#Specify variables for model input
N <- nrow(data)
s <- length(unique(data$Plot))
t <- length(unique(data$Year))
site <- data$Plot
year <- data$Year
y <- data$Germinant_Count

#set directory to save output to
save.dir<-''

#Iterate through each model

for(i in 1:6){
  #covariate matrix
  X <- x[ , names(mods.to.test)[mods.to.test[i,]==1]]
  #fit model in stan
  mod <- stan(file=model_file, data=(list(y=y, N=N, s=s, t=t, site=site, year=year, X=X,
              K=ncol(X))), iter=20000, chains=2, warmup=10000,
              control=list(adapt_delta=0.99,max_treedepth=15),cores=16)
  
  #save out workspace for each model version
  save.image(paste(save.dir, species, '_', model_names[i], '.RData', sep=''))

}
