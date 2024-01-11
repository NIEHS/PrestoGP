set.seed(1234)

#
#
# Multivariate, space-time simulation
# Y(s,t) = X(s) + nu(s,t) + epsilon(s,t)
# Y = Y(i), i = 1,..., k outcomes
#  where:
# X(s) is the spatial only trend,
# nu(s,t) is the spatiotemporally correlated trend, and
# epsilon(s,t) is the spatiotemporal nugget

# To keep the simulation small and tractable to quickly run, we
# simulate 25 spatial locations and 5 temporal locations

n.spatial=25
n.spatial.xy = sqrt(n.spatial)
n.temporal = 5

# Mean trend, X(s): spatiotemporal data, but the trend is spatial only.
# Each outcome has a different combination of true betas
# p = number of potential covariates
p=10

p.nz=4
beta1= c(rep(0,p-p.nz),rep(1,p.nz))
p.nz=2
beta2= c(rep(1,p.nz),rep(0,p-p.nz))
p.nz=3
beta3= c(rep(1,p.nz),rep(0,p-p.nz))


# Combine all of the betas into 1 vector
beta.all <- c(beta1,beta2,beta3)


# correlated predictors
Sigma.X=exp(-rdist(sample(1:p))/3)
X=mvrnorm(n.spatial,rep(0,p),Sigma.X)
# replicate for space-time for a single-outcomes
X.st <- rbind(X,X,X,X,X)
# Create the multivariate X-matrix, for multiple outcomes
# The rows are by outcome (i), space (j), and then time (k)
X.all <- psych::superMatrix(list(X.st,X.st,X.st))


# S-T design matrix with X's repeated is needed

# Calculate the multivariate  mean trend
# The mean trend is spatial only, but now formatted for space-time (i.e. spatial
# covariates are repeated across time)
mean.trend.st <- X.all %*% beta.all


#### Space-Time correlated error ####
# We simulate using a valid parsimonious Matern spatiotemporal cross-covariance function where
# each marginal parameter (spatial range, temporal range, smoothness, marginal variance)
# are all different. The cross-covariance parameters are the average of the respective
# marginal parameters.
# The simulation is done with unconditional Cholesky decomposition


# rho is the correlation between outcome S-T errors, 1 for with itself
rho.vec <- c(0.8, 0.3, 0.5)
rho <- matrix(0,nrow=3,ncol=3)
rho[upper.tri(rho)] <- rho.vec
rho <- rho+t(rho)+diag(1,3)

# nu,  marginal smoothness:
marg.smoothness <- c(0.5,1,1.5)

# Non-spatial error/noise/nugget
# epsilon(s,t) will be simulated based on the nugget parameters,
nuggets <- c(0.5, 0.7, 2)

# marginal variances of the Matern
# The marginal variances are scaled linearly with the nuggets
x.variance <- runif(3,1.5,4)
marg.var <- nuggets + x.variance

# ranges, the true ranges, we will scale coordinates by this, but then use a
# vector of ones in the model/simulation

# outcome1 space and time, outcome 2 space and time, outcome 3 space and time
# Spatial domain is on the unit [0,1]
# Temporal domain is a year [0,365]
ranges <- c(0.8,0.1,0.5)


# set up the true spatial and temporal dimensions
space.dim1 <- seq(0,1,length.out = n.spatial.xy)+rnorm(n.spatial.xy,0,0.001)
space.dim2 <- seq(0,1,length.out = n.spatial.xy)+rnorm(n.spatial.xy,0,0.001)
space.dim3 <- seq(0,1,length.out = n.spatial.xy)+rnorm(n.spatial.xy,0,0.001)
time.dim1 <- seq(0,3,length.out = n.temporal)+rnorm(n.temporal,0,0.0001)
time.dim2 <- seq(0,3,length.out = n.temporal)+rnorm(n.temporal,0,0.0001)
time.dim3 <- seq(0,3,length.out = n.temporal)+rnorm(n.temporal,0,0.0001)

# scale the space and time dimensions according to the true covariance ranges

d11 <- fields::rdist(expand.grid(space.dim1, space.dim1, time.dim1))
d22 <- fields::rdist(expand.grid(space.dim2, space.dim2, time.dim2))
d33 <- fields::rdist(expand.grid(space.dim3, space.dim3, time.dim3))

d12 <- fields::rdist(expand.grid(space.dim1, space.dim1, time.dim1),
                     expand.grid(space.dim2, space.dim2, time.dim2))
d21 <- t(d12)

d13 <- fields::rdist(expand.grid(space.dim1, space.dim1, time.dim1),
                     expand.grid(space.dim3, space.dim3, time.dim3))
d31 <- t(d13)

d23 <- fields::rdist(expand.grid(space.dim2, space.dim2, time.dim2),
                     expand.grid(space.dim3, space.dim3, time.dim3))
d32 <- t(d23)


## Create the correlation matrices
Sigma11 <- marg.var[1] * fields::Matern(d11,range = ranges[1], smoothness = marg.smoothness[1])
Sigma22 <- marg.var[2] * fields::Matern(d22,range = ranges[2], smoothness = marg.smoothness[2])
Sigma33 <- marg.var[3] * fields::Matern(d33,range = ranges[3], smoothness = marg.smoothness[3])

vii <- marg.smoothness[1]
vjj <- marg.smoothness[2]
vij <- (vii+vjj)/2
aii <- 1/ranges[1]
ajj <- 1/ranges[2]
aij <- sqrt((aii^2+ajj^2)/2)
Sigma12 <- rho[1,2] * sqrt(marg.var[1]) * sqrt(marg.var[2]) * aii^vii *
    ajj^vjj * gamma(vij) / (aij^(2*vij) * sqrt(gamma(vii) * gamma(vjj))) *
    fields::Matern(d12, smoothness=vij, alpha=aij)
Sigma21 <- t(Sigma12)

vii <- marg.smoothness[1]
vjj <- marg.smoothness[3]
vij <- (vii+vjj)/2
aii <- 1/ranges[1]
ajj <- 1/ranges[3]
aij <- sqrt((aii^2+ajj^2)/2)
Sigma13 <- rho[1,3] * sqrt(marg.var[1]) * sqrt(marg.var[3]) * aii^vii *
    ajj^vjj * gamma(vij) / (aij^(2*vij) * sqrt(gamma(vii) * gamma(vjj))) *
    fields::Matern(d13, smoothness=vij, alpha=aij)
Sigma31 <- t(Sigma13)

vii <- marg.smoothness[2]
vjj <- marg.smoothness[3]
vij <- (vii+vjj)/2
aii <- 1/ranges[2]
ajj <- 1/ranges[3]
aij <- sqrt((aii^2+ajj^2)/2)
Sigma23 <- rho[2,3] * sqrt(marg.var[2]) * sqrt(marg.var[3]) * aii^vii *
    ajj^vjj * gamma(vij) / (aij^(2*vij) * sqrt(gamma(vii) * gamma(vjj))) *
    fields::Matern(d23, smoothness=vij, alpha=aij)
Sigma32 <- t(Sigma23)

#Combine into the super Multivariate covariance matrix
Sigma.All <- rbind( cbind(Sigma11,Sigma12,Sigma13),
                    cbind(Sigma21,Sigma22,Sigma23),
                    cbind(Sigma31,Sigma32,Sigma33)
  )

# Cholesky decomposition
L.C <- chol(Sigma.All)
# transpose
L.Ct <- t(L.C)

locs.list <- list()
locs.list[[1]] <- as.matrix(expand.grid(space.dim1, space.dim1, time.dim1))
locs.list[[2]] <- as.matrix(expand.grid(space.dim2, space.dim2, time.dim2))
locs.list[[3]] <- as.matrix(expand.grid(space.dim3, space.dim3, time.dim3))
