library(PrestoGP)

set.seed(7919)
load("tests/testthat/multivariate_sim.Rdata")
model <- new("MultivariateSpatialModel")
model <- prestogp_fit(model, Y_train, X_train, locs_train, verbose = TRUE)
prediction <- prestogp_predict(model, X_test, locs_test)
means <- prediction[[1]]
mse <- mean((Y_test-means)^2)
mse
sds <- prediction[[2]]
#simple_model = cv.ncvreg(as.matrix(X_train), Y_train, family = "gaussian",
#                                        penalty = "SCAD",dfmax=100,returnX = FALSE)
#simple_model_means <- predict(simple_model, X = X_test)
#simple_mse <- mean((Y_test-simple_model_means)^2)