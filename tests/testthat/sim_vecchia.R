set.seed(1212)

p <- 10 # number of predictors for each response
p.nz <- 4 # number of nonzero predictors for each y
n.spatial.xy <- 50 # number of spatial coordinates per dimension

library(MASS)
library(fields)

beta1 <- c(rep(1, p.nz), rep(0, p - p.nz))

Sigma.X <- exp(-rdist(sample(1:p)) / 3)
X <- mvrnorm(n.spatial.xy^2, rep(0, p), Sigma.X)
mean.trend.st <- as.vector(X %*% beta1)

marg.smoothness <- 0.5 + rnorm(1, 0.1, 0.05)
nuggets <- runif(1, 0.5, 2)
x.variance <- runif(1, 1.5, 4)
marg.var <- nuggets + x.variance
ranges <- runif(1, 0.5, 1.2)

params.all <- c(x.variance, ranges, marg.smoothness, nuggets)


loc1 <- seq(0, 1, length.out = n.spatial.xy) + rnorm(n.spatial.xy, 0, 0.001)
loc2 <- seq(0, 1, length.out = n.spatial.xy) + rnorm(n.spatial.xy, 0, 0.001)
locs <- as.matrix(expand.grid(loc1, loc2))

Sigma.All <- marg.var * Matern(rdist(locs),
   range = ranges,
   smoothness = marg.smoothness
)
L.C <- chol(Sigma.All)

st.error <- as.vector(rnorm(n.spatial.xy^2) %*% L.C)
nug.error <- nuggets * rnorm(n.spatial.xy^2)
y <- mean.trend.st + st.error + nug.error

rm(
   p, p.nz, n.spatial.xy, beta1, Sigma.X, mean.trend.st, loc1, loc2, L.C,
   st.error, nug.error, ranges, Sigma.All, nuggets, marg.smoothness, marg.var,
   x.variance
)
