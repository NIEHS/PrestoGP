set.seed(1212)

p <- 10 # number of predictors for each response
p.nz <- 4 # number of nonzero predictors for each y
n.spatial.xy <- 15 # number of spatial coordinates per dimension

beta1 <- c(rep(1, p.nz), rep(0, p - p.nz))

Sigma.X <- exp(-rdist(sample(1:p)) / 3)
X <- mvrnorm(n.spatial.xy^3, rep(0, p), Sigma.X)
mean.trend.st <- as.vector(X %*% beta1)

marg.smoothness <- 0.5 + rnorm(1, 0.1, 0.05)
nuggets <- runif(1, 0.5, 2)
x.variance <- runif(1, 1.5, 4)
marg.var <- nuggets + x.variance
ranges1 <- runif(1, 0.5, 1.2)
ranges2 <- runif(1, 5, 10)

params.all <- c(marg.var, ranges1, ranges2, marg.smoothness, nuggets)


loc1 <- seq(0, 1, length.out = n.spatial.xy) + rnorm(n.spatial.xy, 0, 0.001)
loc2 <- seq(0, 1, length.out = n.spatial.xy) + rnorm(n.spatial.xy, 0, 0.001)
loc3 <- seq(0, 1, length.out = n.spatial.xy) + rnorm(n.spatial.xy, 0, 0.001)
locs <- as.matrix(expand.grid(loc1, loc2, loc3))

locss <- locs
locss[, 1:2] <- locss[, 1:2] / ranges1
locss[, 3] <- locs[, 3] / ranges2

Sigma.All <- marg.var * Matern(rdist(locss), range = 1,
  smoothness = marg.smoothness)
L.C <- chol(Sigma.All)

st.error <- as.vector(rnorm(n.spatial.xy^3) %*% L.C)
nug.error <- sqrt(nuggets) * rnorm(n.spatial.xy^3)
y <- mean.trend.st + st.error + nug.error
