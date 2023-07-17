# LURK-Vecchia

Simultaneous variable seletion and estimation of LUR models with spatiotemporally correlated errors that is scalable for big data

# Citation
Messier, K.P. & Katzfuss, M. (2020). Scalable penalized spatiotemporal land-use regression for ground-level nitrogen dioxide. https://arxiv.org/abs/2005.09210

# Installation 
Based on the package GPvecchia (https://cran.r-project.org/web/packages/GPvecchia/index.html) and
ncvreg (https://cran.r-project.org/web/packages/ncvreg/index.html). 
Install the following R packages: GPvecchia,ncvreg, fields, Matrix

Download LURK_Functions.R to your local machine and use R to source:
```
source("/Location/of/your/functions/LURK_Functions.R") 
``` 

# Data

## NO<sub>2</sub>
The NO<sub>2</sub> analysis data are stored in a single csv file in the [data](https://github.com/NIEHS/LURK-Vecchia/tree/master/data) subfolder. 

```
US_ST_NO2_Data.csv
```

Columns 1 - 6 are latitude, longitude, date, NO2 concentration, and time in fractional years, respectively. Columns 7 - 145 make up the design matrix of covariates. The covariate names are followed by the distance hyperparameter (i.e. buffer, decay range).
To see all of the field names, call: 
```
colnames(US_ST_NO2_Data)
```

## Simulations

The simulation analysis data are stored in Excel files in the [data](https://github.com/NIEHS/LURK-Vecchia/tree/master/data) subfolder. 
There is an Excel file for each spatiotemporal parameter that is varied: nugget-to-sill ratio, total variance (sill), spatial range, and temporal range.
For each spatiotemporal parameter there is a training set and test set, stored in separate files. The training and test sets are based on the same
simulation.

```
Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx
Vecchia_ST_Simulation_NUG2SILL5_test_20191107.xlsx
Vecchia_ST_Simulation_VARIANCE5_20191107.xlsx
Vecchia_ST_Simulation_VARIANCE5_test_20191107.xlsx
Vecchia_ST_Simulation_SPATIAL_RANGES5_20191107.xlsx
Vecchia_ST_Simulation_SPATIAL_RANGES5_test_20191107.xlsx
Vecchia_ST_Simulation_TEMPORAL_RANGES5_20191107.xlsx
Vecchia_ST_Simulation_TEMPORAL_RANGES5_test_20191107.xlsx
```
Each Excel file contains a sheet for each simulation scenario - that is there are 20 sheets for 20 variations of the spatiotemporal parameter. 
Each sheet is simply labeled "1", "2",..., etc. for the simulation scenario. 
The first sheet also contains data for the true observations without measurement error and the covariates. 
Hence, in each simulation you see differences for the first sheet. 

Here we go through an **example** of reading in data and assigning variables to help demonstrate the data format:
```
data1 <- as.data.frame(read_excel("Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx", 
                                  sheet = "1", col_names = FALSE))
```
The true data, without measurement error is read in:

```
# y true data
Y.true <- as.matrix(data1[,1])
```
Columns 2 through 21 on the first sheet are the 20 simulation scenarios of y-observed (i.e. with measurement error). 
The covariates are read in.
```
# covariates
X <- as.matrix(data1[,22:144])
```
There is a separate sheet for the spatial coordiantes:
```
# longitude and latitude
xyt <- as.matrix(read_excel("Vecchia_ST_Simulation_NUG2SILL5_20191107.xlsx", 
                            sheet = "xy", col_names = FALSE))

```
The same information is also read in for the test set, which are 1000 other space-time observations from the same
respective simulation. 
```
xy_test <- as.matrix(read_excel("Vecchia_ST_Simulation_NUG2SILL5_test_20191107.xlsx", 
                                sheet = "1", col_names = FALSE))

test_xyt <- as.matrix(read_excel("Vecchia_ST_Simulation_NUG2SILL5_test_20191107.xlsx", 
                                 sheet = "xy", col_names = FALSE))
y_test <- as.matrix(xy_test[,1])

X_test <- as.matrix(xy_test[,22:144])
```

The main simulations go through nested loops: The outer loop (i) goes through each parameter variation (e.g. differing nugget-to-sill ratio). The inner loop (j)
goes through random simulations of the same parameter values. Example below:
```
for (i in 1:20){
  
  ### Read in the simulated data
  Ydata <- read_excel("Vecchia_ST_Simulation_TEMPORAL_RANGES5_20191107.xlsx", 
                      sheet = as.character(i), col_names = FALSE)
  
  test_data <- read_excel("Vecchia_ST_Simulation_TEMPORAL_RANGES5_test_20191107.xlsx", 
                          sheet = as.character(i), col_names = FALSE)
```
The following if-else statement deals with the true data in column 1 of sheet 1.

```
  
  if (i == 1){
    Y.all <- Ydata[,2:21]
    test_y <- test_data[,2:21]
  }else{
    Y.all <- Ydata[,1:20]
    test_y <- test_data[1:20]
  }
  
````
Then the inner loop (j) goes through random realizations
```  
  for (j in 1:P){
    

    ########## Observed data for simulation j ############
    Y.obs <- data.matrix(Y.all[,j])
    n=length(Y.obs)
    test_y_j <- data.matrix(test_y[,j])
```


# Running the Simulation Analysis 

## A runnable simulation example
The beginning of the simulation contains one iteration of the paper's full simulations. This is
included so that a user can run a simulation and produce the results in a timely manner. Due to the 
number of simulations, the full paper simulations will take a long time. Note that the main simulations were sub-divided
and ran on high-performance computing clusters.

The runnable simulation is the first section of Analysis_Simulation.R. The header of the section starts with:
```
####################################################################################################
### SECTION 1
### A short simulation scenario that can be run in a few minutes
```

## Replicating the paper simulations

The paper simulations have a section for each scenario. 
- Section 2: Variations of the temporal range
- Section 3: Variations of the spatial range
- Section 4: Variations of the nugget-to-sill ratio
- Section 5: Variations of the total variance 
- Section 6: Convergence of estimated LUR coefficients (betas) with increasing m in the Vecchia approximation.

# Running the NO<sub>2</sub> Analysis 

With the LURK_Functions sourced and the data accessible, run Analysis_NO2.R to replicate the 10-fold cross-validation
NO<sub>2</sub> from the paper.
