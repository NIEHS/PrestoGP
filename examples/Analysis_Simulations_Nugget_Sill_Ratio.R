library(LURKVecchia)

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

####################################################################################################
####################################################################################################
### SECTION 4
### Nugget to Sill Ratio simulations (20 simulations for 20 variations of the nug-to-sill ratio ) ###########
####################################################################################################
####################################################################################################

###### Load the simulation data
data1 <- as.data.frame(read_excel("Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx", 
                                  sheet = "1", col_names = FALSE))

# y true data
Y.true <- as.matrix(data1[,1])
# covariates
X <- as.matrix(data1[,22:144])

# longitude and latitude
xyt <- as.matrix(read_excel("Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx", 
                            sheet = "xy", col_names = FALSE))


# Read in the test set data 
xy_test <- as.matrix(read_excel("Vecchia_ST_Simulation_NUG2SILL5_test_20191107.xlsx", 
                                sheet = "1", col_names = FALSE))

test_xyt <- as.matrix(read_excel("Vecchia_ST_Simulation_NUG2SILL5_test_20191107.xlsx", 
                                 sheet = "xy", col_names = FALSE))
y_test <- as.matrix(xy_test[,1])

X_test <- as.matrix(xy_test[,22:144])



### Pre-allocate the validation statistic vectors 
P = 20
Local.Kriging.crps = matrix(NA,nrow= P * 20)
Local.Kriging.logs = matrix(NA,nrow= P * 20)
Local.Kriging.MSE = matrix(NA,nrow= P * 20)
LUR.iid.MSE = matrix(NA,nrow= P * 20)
LUR.iid.crps = matrix(NA,nrow= P * 20)
LUR.iid.logs = matrix(NA,nrow= P * 20)
LURK.Local.crps = matrix(NA,nrow= P * 20)
LURK.Local.logs = matrix(NA,nrow= P * 20)
LURK.Local.MSE = matrix(NA,nrow= P * 20)
LURK.Full.MSE = matrix(NA,nrow= P * 20)
LURK.Full.logs = matrix(NA,nrow= P * 20)
LURK.Full.crps = matrix(NA,nrow= P * 20)
LURK.Vecchia.MSE = matrix(NA,nrow= P * 20)
LURK.Vecchia.logs = matrix(NA,nrow= P * 20)
LURK.Vecchia.crps = matrix(NA,nrow= P * 20)

LUR.iid.beta.estimates = matrix(NA,nrow = ncol(X)+1,ncol = P * 20)
LURK.Vecchia.beta.estimates = matrix(NA,nrow = ncol(X)+1,ncol = P * 20)


########################
########################
### Begin Main Loop #####

iter=1
for (i in 1:20){
  
  ### Read in the simulated data
  Ydata <- read_excel("Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx", 
                      sheet = as.character(i), col_names = FALSE)
  
  test_data <- read_excel("Vecchia_ST_Simulation_NUG2SILL5_test_20191107.xlsx", 
                          sheet = as.character(i), col_names = FALSE)
  
  if (i == 1){
    Y.all <- Ydata[,2:21]
    test_y <- test_data[,2:21]
  }else{
    Y.all <- Ydata[,1:20]
    test_y <- test_data[1:20]
  }
  
  
  for (j in 1:P){
    
    
    ########## Observed data for simulation j ############
    Y.obs <- data.matrix(Y.all[,j])
    n=length(Y.obs)
    test_y_j <- data.matrix(test_y[,j])
    
    
    ########## Local-Kriging ###############
    print(i)
    print("Local-Kriging")
    
    Local.Kriging.cov.param = ST_Krig_Param_Avg(Y.obs,xyt,p = 116,k = 10)
    
    
    Local.Kriging.pred <- Kr_pred(test_xyt,xyt,Y.obs,Local.Kriging.cov.param,25)  
    Local.Kriging.mean<- Local.Kriging.pred$Kr.prediction 
    Local.Kriging.sd <- sqrt(Local.Kriging.pred$Kr.Var)
    
    
    Local.Kriging.MSE[iter] <- mean((Local.Kriging.mean- test_y_j)^2)
    Local.Kriging.logs[iter] <- mean(logs_norm(test_y_j,Local.Kriging.mean,Local.Kriging.sd))    
    Local.Kriging.crps[iter] <- mean(crps_norm(test_y_j,Local.Kriging.mean,Local.Kriging.sd))    
    
    ########## LUR-iid (LUR.iid) ###############
    
    print(i)
    print("LUR.iid")
    LUR.iid.fit=cv.ncvreg(as.matrix(X),Y.obs,family = "gaussian",
                          penalty = "SCAD",returnX = FALSE,dfmax = 100)
    
    idmin <- which(LUR.iid.fit$lambda == LUR.iid.fit$lambda.min)
    semin <- LUR.iid.fit$cve[idmin] + LUR.iid.fit$cvse[idmin]
    lambda.1se <- max(LUR.iid.fit$lambda[LUR.iid.fit$cve<=semin])
    lambda.1se.idx <- which(LUR.iid.fit$lambda==lambda.1se)
    
    LUR.iid.pred <- predict(LUR.iid.fit,X = as.matrix(X_test),lambda=lambda.1se) 
    LUR.iid.beta.estimates[,iter] <- coef(LUR.iid.fit,lambda=lambda.1se)
    beta.in <- LUR.iid.beta.estimates!=0
    # Get the SCAD prediction variance assuming LR assumptions  
    tau2.LUR.iid <- (1/(n-sum(LUR.iid.beta.estimates[,iter]!=0))) * norm(y_test - LUR.iid.pred,type = "F")^2
    
    LUR.iid.MSE[iter] = mean((LUR.iid.pred-test_y_j)^2)
    LUR.iid.crps[iter] <- mean(crps_norm(test_y_j,LUR.iid.pred, rep(tau2.LUR.iid,length(LUR.iid.pred))))
    LUR.iid.logs[iter] <- mean(logs_norm(test_y_j,LUR.iid.pred, rep(tau2.LUR.iid,length(LUR.iid.pred))))
    
    ########## LURK-Local ###### ###############
    
    # de-trend
    print(i)
    print("LURK.Local")
    LUR.iid.train.pred <- predict(LUR.iid.fit,X = as.matrix(X),lambda=lambda.1se)
    res.train <- Y.obs- LUR.iid.train.pred
    
    
    LURK.Local.cov.param = ST_Krig_Param_Avg(res.train,xyt,p = 116,k = 10)
    
    
    LURK.Local.pred <- Kr_pred(test_xyt,xyt,res.train,LURK.Local.cov.param,25)  
    LURK.Local.mean<- LURK.Local.pred$Kr.prediction + LUR.iid.pred
    LURK.Local.sd <- sqrt(LURK.Local.pred$Kr.Var)
    
    
    
    LURK.Local.MSE[iter] <- mean((LURK.Local.mean- test_y_j)^2)
    LURK.Local.logs[iter] <- mean(logs_norm(test_y_j,LURK.Local.mean,LURK.Local.sd))    
    LURK.Local.crps[iter] <- mean(crps_norm(test_y_j,LURK.Local.mean,LURK.Local.sd))    
    
    
    ########## LURK-Full ###### ###############
    
    # LURK-Full estimation of the coviarance and beta parameters
    results.LURK.Full <- LURK_Full(Y = Y.obs,X = as.matrix(X), locs = xyt)
    
    Full.beta.estimates <- results.LURK.Full$beta  
    Full.theta <- results.LURK.Full$covparam
    
    # Full prediction at new locations
    full.trend.prediction <- predict(results.LURK.Full$Full.SCAD.fit,X = X_test, which = results.LURK.Full$lambda.1se.idx) 
    # Full trend prediction at observed data
    LURK.Full.hat <- predict(results.LURK.Full$Full.SCAD.fit,X = X, which = results.LURK.Full$lambda.1se.idx)
    
    ### Test set prediction
    # calculate observed residuals
    res=Y.obs - LURK.Full.hat
    
    # get the scaled S-T coordinates
    locs.train.scaled = cbind(xyt[,1]/Full.theta[2], 
                              xyt[,2]/Full.theta[2], xyt[,3]/Full.theta[3])  
    
    locs.test.scaled = cbind(test_xyt[,1]/Full.theta[2], 
                             test_xyt[,2]/Full.theta[2], test_xyt[,3]/Full.theta[3])  
    
    # calculate distance matrices
    dist.oo <- fields::rdist(locs.train.scaled,locs.train.scaled)
    dist.op <- fields::rdist(locs.train.scaled,locs.test.scaled)
    
    # calculate the covariance matrices 
    Sigma.oo <- Full.theta[4]*diag(n) + 
      Full.theta[1]*fields::Exponential(dist.oo,range=1)
    Sigma.op <- Full.theta[1] * fields::Exponential(dist.op,range = 1)
    
    # mean of residuals
    mean_trend <- mean(res)
    
    # Kriging prediction
    Kr.prediction <- mean_trend + t(Sigma.op) %*% solve(Sigma.oo,res - mean_trend)
    
    # Kriging standard deviation
    Full.sd = sqrt(diag(Full.theta[1]-t(Sigma.op) %*% solve(Sigma.oo, Sigma.op)))
    
    # Add 
    Full.mean<- Kr.prediction + full.trend.prediction
    
    
    
    LURK.Full.MSE[iter] <- mean((Full.mean - test_y_j)^2)
    LURK.Full.logs[iter] <- mean(logs_norm(test_y_j,Full.mean,Full.sd))    
    LURK.Full.crps[iter] <- mean(crps_norm(test_y_j,Full.mean,Full.sd))   
    
    
    ########## LURK-Vecchia ###### ###############
    m=25
    
    # LURK-Vecchia estimation of the coviarance and beta parameters
    results.LURK.Vecchia <- LURK_Vecchia(Y = Y.obs,X = X,locs = xyt, m = m)
    Vecchia.beta.estimates <- results.LURK.Vecchia$beta  
    LURK.Vecchia.beta.estimates[,i] <- Vecchia.beta.estimates
    Vecchia.theta <- results.LURK.Vecchia$covparam
    
    # Vecchia prediction at new locations
    LURK.Vecchia.Trend.Pred <- predict(results.LURK.Vecchia$Vecchia.SCAD.fit,X = X_test, which = results.LURK.Vecchia$lambda.1se.idx) 
    # Vecchia trend prediction at observed data
    LURK.Vecchia.hat <- predict(results.LURK.Vecchia$Vecchia.SCAD.fit,X = X, which = results.LURK.Vecchia$lambda.1se.idx)
    
    # Test set prediction
    res=Y.obs - LURK.Vecchia.hat
    
    locs.train.scaled = cbind(xyt[,1]/Vecchia.theta[2], 
                              xyt[,2]/Vecchia.theta[2], xyt[,3]/Vecchia.theta[3])  
    
    locs.test.scaled = cbind(test_xyt[,1]/Vecchia.theta[2], 
                             test_xyt[,2]/Vecchia.theta[2], test_xyt[,3]/Vecchia.theta[3])  
    
    
    vec.approx.test=vecchia_specify(locs.train.scaled,m,locs.pred=locs.test.scaled)
    
    
    ## carry out prediction
    pred=vecchia_prediction(res,vec.approx.test,c(Vecchia.theta[1],1,0.5),Vecchia.theta[4])
    
    # returns a list with elements mu.pred,mu.obs,var.pred,var.obs,V.ord
    Vec.mean=pred$mu.pred+LURK.Vecchia.Trend.Pred
    Vec.sds=sqrt(pred$var.pred)
    
    LURK.Vecchia.MSE[iter] <- mean((Vec.mean - test_y_j)^2)
    LURK.Vecchia.logs[iter] <- mean(logs_norm(test_y_j,Vec.mean,Vec.sds))    
    LURK.Vecchia.crps[iter] <- mean(crps_norm(test_y_j,Vec.mean,Vec.sds))   
    
    ################################################## end of iteration
    
    
    
    iter = iter+1
    
  } # P 
  
  save(Local.Kriging.MSE,Local.Kriging.crps,Local.Kriging.logs,
       LUR.iid.MSE,LUR.iid.crps,LUR.iid.logs,LUR.iid.beta.estimates,
       LURK.Local.MSE,LURK.Local.crps,LURK.Local.logs,
       LURK.Full.MSE,LURK.Full.crps,LURK.Full.crps,
       LURK.Vecchia.beta.estimates,LURK.Vecchia.crps,LURK.Vecchia.logs,LURK.Vecchia.beta.estimates,
       file = "Simulation_Results_by_NUG2SILL.RData")
  
  
  
}#Outer loop (Excel Sheets) 
#### END OF THE SECTION FOR THE SIMULATIONS BY NUGGET-TO-SILL RATIO