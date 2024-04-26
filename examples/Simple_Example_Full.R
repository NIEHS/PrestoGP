library(PrestoGP)

# For reading in Excel data 
library(readxl)

###### Load the simulation data
data1 <- as.data.frame(read_excel("data/Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx", 
                                  sheet = "1", col_names = FALSE))

# y true data
Y.true <- as.matrix(data1[,1])
# covariates
X <- as.matrix(data1[,22:144])

# longitude and latitude
xyt <- as.matrix(read_excel("data/Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx", 
                            sheet = "xy", col_names = FALSE))


# Read in the test set data 
xy_test <- as.matrix(read_excel("data/Vecchia_ST_Simulation_NUG2SILL5_test_20191107.xlsx", 
                                sheet = "1", col_names = FALSE))

test_xyt <- as.matrix(read_excel("data/Vecchia_ST_Simulation_NUG2SILL5_test_20191107.xlsx", 
                                 sheet = "xy", col_names = FALSE))
y_test <- as.matrix(xy_test[,1])

X_test <- as.matrix(xy_test[,22:144])



####################
sim.iter = 3

### Read in the simulated data
Ydata <- read_excel("data/Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx", 
                    sheet = as.character(sim.iter), col_names = FALSE)

test_data <- read_excel("data/Vecchia_ST_Simulation_NUG2SILL5_test_20191107.xlsx", 
                        sheet = as.character(sim.iter), col_names = FALSE)

Y.all <- Ydata[,1:20]
test_y <- test_data[1:20]


########## Observed data for simulation j ############
Y.obs <- data.matrix(Y.all[,1])
n=length(Y.obs)
test_y_j <- data.matrix(test_y[,1])


########## LURK-Vecchia ###### ###############

full_model <- new("SpatiotemporalFullModel")
full_model <- lurk_fit(full_model, Y.obs, X, xyt, verbose = TRUE)

prediction <- lurk_predict(full_model, X_test, test_xyt)
Vec.mean <- prediction[[1]]
Vec.sds <- prediction[[2]]
LURK.Full.beta.estimates <- full_model@beta
Full.theta <- full_model@covparams

LURK.Full.MSE <- mean((Vec.mean - test_y_j)^2)
