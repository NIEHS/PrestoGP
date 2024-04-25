library(PrestoGP)

set.seed(7919)
load("tests/testthat/sim_spatial.Rdata")
model <- new("SpatialModel", n_neighbors=4)
model <- prestogp_fit(model, Y_train, X_train, locs_train, verbose = TRUE)
prediction <- prestogp_predict(model, X_test, locs_test, m = 4)
means <- prediction[[1]]
#mse <- mean((Y_test-means)^2)
mse <- crossprod((Y_test-means)) / (nrow(Y_test)-colSums(model@beta))
diag(mse)
plot(means, Y_test)
msds <- mean(prediction[[2]])
msds