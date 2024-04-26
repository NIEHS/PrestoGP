library(PrestoGP)

set.seed(7919)
load("tests/testthat/sim_spatial.Rdata")
model <- new("MultivariateSpatialModel")
model <- prestogp_fit(model, Y_train, X_train, locs_train, verbose = TRUE)
prediction <- prestogp_predict(model, X_test, locs_test)
means <- prediction[[1]]
mse <- mean((Y_test-means)^2)
mse
sds <- mean(prediction[[2]])
sds
#simple_model = cv.ncvreg(as.matrix(X_train), Y_train, family = "gaussian",
#                                        penalty = "SCAD",dfmax=100,returnX = FALSE)
#simple_model_means <- predict(simple_model, X = X_test)
#simple_mse <- mean((Y_test-simple_model_means)^2)

plot(means, Y_test)
plot(means, prediction[[2]])

df.plot <- data.frame("x.coord" = round(locs_test[,1], 0),
                      "y.coord" = round(locs_test[,2], 0),
                      "time" = as.factor(rep(0, nrow(locs_test))),
                      #"time" = as.factor(xyt.grid[,3]),
                      "var1" = Y_test,
                      "var1.pred" = means,
                      "sds" = prediction[[2]])
#Var 1
ggplot(df.plot,aes(x.coord,y.coord,fill = var1))+geom_raster(interpolate = T)+
  facet_wrap(.~time)+scale_fill_viridis_c(direction = -1, option = "A")

ggplot(df.plot,aes(x.coord,y.coord,fill = var1.pred))+geom_raster(interpolate = T)+
  facet_wrap(.~time)+scale_fill_viridis_c(direction = -1, option = "A")

ggplot(df.plot,aes(x.coord,y.coord,fill = sds))+geom_raster(interpolate = T)+
  facet_wrap(.~time)+scale_fill_viridis_c(direction = -1, option = "A")