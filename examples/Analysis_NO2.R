
###### Libraries needed in Vecchia ###########
library(GPvecchia)
library(Matrix)
library(fields)
### SCAD
library(ncvreg)

# For reading in Excel data 
library(readxl)
# Proper scoring rules (log-score and CRPS)
library(scoringRules)


# Example path - needs to be changed for your local machine
#source("R/LURK_Functions.R")

###### Set the Working directory, if needed
#setwd()


### DATA #####
US_NO2_Data<- read.csv("data/US_ST_NO2_Data.csv")

lat       <- US_NO2_Data$Latitude # Latitude
lon       <- US_NO2_Data$Longitude # Longitude
logNO2     <-  log(US_NO2_Data$Y) # ozone data
time      <- US_NO2_Data$YearFrac    # time in fractional years
N         <- length(logNO2)
locs = cbind(lon,lat,time) # Coordinates in R3 (x,y,t)

X.Design <- US_NO2_Data[,7:145] #Design matrix
X.Design <- scale(X.Design)
ID <- US_NO2_Data$ID
# Set the seed and create foldids - seed is necessary for repeatability
set.seed(1)
foldid <- sample(1:10,N,replace=TRUE)

Local.Kriging.MSE <- matrix(nrow = 10)
Local.Kriging.crps <- matrix(nrow = 10)
Local.Kriging.logs <- matrix(nrow = 10)

LUR.iid.MSE <- matrix(nrow = 10)
LUR.iid.crps <- matrix(nrow = 10)
LUR.iid.logs <- matrix(nrow = 10)

LURK.Local.MSE <- matrix(nrow = 10)
LURK.Local.crps <- matrix(nrow = 10)
LURK.Local.logs <- matrix(nrow = 10)


LURK.Vecchia.MSE <- matrix(nrow = 10)
LURK.Vecchia.crps <- matrix(nrow = 10)
LURK.Vecchia.logs <- matrix(nrow = 10)


########################
### Begin Main Loop - 10 fold cross-validation of NO2 dataset #####
for (i in 1:10){
  
  
  # Set up the training and test data
  idx.train <- foldid != i
  idx.test  <- foldid == i
  
  #### FOR TESTING PURPOSE - use test set size
  # idx.train <- idx.test
  ###
  locs.train <- locs[idx.train,]
  locs.test  <- locs[idx.test,]
  X.train <- X.Design[idx.train,]
  X.test  <- X.Design[idx.test,]
  Y.train <- logNO2[idx.train]
  Y.test  <- logNO2[idx.test]
  n=sum(idx.train)
  
  ################### 
  ### Pure Geostatistical estimation (Local Simple Kriging)
  ###################
  print(i)
  print("OK")
  
  
  OK.cov.param = ST_Krig_Param_Avg(Y.train,locs.train,p = 500,k = 10)
  
  #(new_coords,obs_coords,Y_obs,cov.pars,NN)
  OK.pred <- Kr_pred(locs.test,locs.train,Y.train,OK.cov.param,25)  
  OK.mean<- OK.pred$Kr.prediction
  OK.sd <- sqrt(OK.pred$Kr.Var + OK.cov.param[4])
  
  Local.Kriging.MSE[i] <- mean((OK.mean- Y.test)^2)
  Local.Kriging.logs[i] <- mean(logs_norm(Y.test,OK.mean,OK.sd))    
  Local.Kriging.crps[i] <- mean(crps_norm(Y.test,OK.mean,OK.sd))    
  
  ######## iid SCAD #############
  print(i)
  print("SCAD")
  SCAD.fit=cv.ncvreg(as.matrix(X.train),Y.train,family = "gaussian",
                     penalty = "SCAD",returnX = FALSE)
  
  idmin <- which(SCAD.fit$lambda == SCAD.fit$lambda.min)
  semin <- SCAD.fit$cve[idmin] + SCAD.fit$cvse[idmin]
  lambda.1se <- max(SCAD.fit$lambda[SCAD.fit$cve<=semin])
  lambda.1se.idx <- which(SCAD.fit$lambda==lambda.1se)
  
  SCAD.pred <- predict(SCAD.fit,X = as.matrix(X.test),lambda=lambda.1se) 
  SCAD.beta.estimates <- coef(SCAD.fit,lambda=lambda.1se)
  beta.in <- SCAD.beta.estimates!=0
  # Get the SCAD prediction variance assuming LR assumptions  
  tau2.SCAD <- (1/(n-sum(SCAD.beta.estimates!=0))) * norm(Y.test - SCAD.pred,type = "F")^2
  
  LUR.iid.MSE[i] = mean((SCAD.pred-Y.test)^2)
  LUR.iid.crps[i] <- mean(crps_norm(Y.test,SCAD.pred, rep(sqrt(tau2.SCAD),length(SCAD.pred))))
  LUR.iid.logs[i] <- mean(logs_norm(Y.test,SCAD.pred, rep(sqrt(tau2.SCAD),length(SCAD.pred))))
  
  ######## Kriging with an external drift  #############
  # de-trend
  print(i)
  print("KED")
  SCAD.train.pred <- predict(SCAD.fit,X = as.matrix(X.train),lambda=lambda.1se)
  res.train <- Y.train- SCAD.train.pred
  
  KED.cov.param = ST_Krig_Param_Avg(res.train,locs.train,p = 500,k = 10)
  
  
  KED.pred <- Kr_pred(locs.test,locs.train,res.train,KED.cov.param,25)  
  KED.mean<- KED.pred$Kr.prediction + SCAD.pred
  KED.sd <- sqrt(KED.pred$Kr.Var + KED.cov.param[4])
  
  
  
  LURK.Local.MSE[i] <- mean((KED.mean- Y.test)^2)
  LURK.Local.logs[i] <- mean(logs_norm(Y.test,KED.mean,KED.sd))    
  LURK.Local.crps[i] <- mean(crps_norm(Y.test,KED.mean,KED.sd))    
  
  
  
  ##############################################
  #########LURK-Vecchia ########################
  ##############################################
  
  model <- SpatiotemporalModel()
  model <- lurk_fit(model, as.matrix(Y.train), X.train, locs.train, verbose = TRUE)
  
  prediction <- lurk_predict(model, X.test, locs.test)
  Vec.mean <- prediction[[1]]
  Vec.sds <- prediction[[1]]
  
  LURK.Vecchia.MSE[i] <- mean((Vec.mean - Y.test)^2)
  #See https://cran.r-project.org/web/packages/scoringRules/vignettes/gettingstarted.html
  LURK.Vecchia.logs[i] <- mean(logs_norm(Y.test,Vec.mean,Vec.sds)) #logarithmic score
  LURK.Vecchia.crps[i] <- mean(crps_norm(Y.test,Vec.mean,Vec.sds)) #continuous ranked probability score 
  
  
  
}#Outer loop (10 fold cross-validation)



save(Local.Kriging.crps,Local.Kriging.logs,Local.Kriging.MSE,
     LUR.iid.crps,LUR.iid.logs,LUR.iid.MSE,
     LURK.Local.crps,LURK.Local.logs,LURK.Local.MSE,
  LURK.Vecchia.crps,LURK.Vecchia.logs,LURK.Vecchia.MSE,
      file = "ValStats_NO2_10fold_Crossvalidation.RData")


