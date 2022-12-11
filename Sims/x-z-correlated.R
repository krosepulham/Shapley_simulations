#setup##########################################################################
library(foreach)
library(doParallel)
library(shapleyAIC)
#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

rm(list = ls())
set.seed(22116)
AICvectors <- list()
rsqvectors <- list()
nsim  <- 1000
nsamp <- 50
beta_true <- c(1,1,1,1,0,0,0,0)
start_time = Sys.time()
output_file <- "Results/x-z-correlated.RData"

#log start of simulation########################################################
sink(file="log.txt",append=TRUE,type="output")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Beginning simulation with:\n")
cat("nsim = ",nsim,"\n")
cat("nsamp = ",nsamp,"\n")
cat("beta_true = ",beta_true,"\n")
cat("...")
sink()

#simulation loop################################################################
#for loop simulates the models and stores the shapley values
start_time = Sys.time()
finalMatrix <- foreach(i=1:nsim, .combine=rbind) %dopar% {
  X1 <- rnorm(nsamp)
  X2 <- rnorm(nsamp)
  X3 <- rnorm(nsamp)
  X4 <- rnorm(nsamp)
  Z1 <- X1 + rnorm(nsamp,mean = 0,sd=1/4)
  Z2 <- X2 + rnorm(nsamp,mean = 0,sd=1/4)
  Z3 <- X3 + rnorm(nsamp,mean = 0,sd=1/4)
  Z4 <- X4 + rnorm(nsamp,mean = 0,sd=1/4)

  #put together into design matrix
  X <- cbind(X1,X2,X3,X4,Z1,Z2,Z3,Z4)

  #calculate Y
  Y <- X%*%cbind(beta_true)+rnorm(nsamp)

  #calculate variables shapley values
  Yshap <- shapleyAIC::shapley(Y,X)
  #AICvectors[[i]]<- Yshap$AIC_shapleys
  #rsqvectors[[i]]<- Yshap$rsq_shapleys

  tempMatrix <- c(Yshap$AIC_shapleys,Yshap$rsq_shapleys)
  tempMatrix #Equivalent to finalMatrix = rbind(finalMatrix, tempMatrix)
}
#old code from before parallelizing
# for(i in 1:nsim){
#   #randomly generate explanatory covariates
#   X1 <- rnorm(nsamp)
#   X2 <- rnorm(nsamp)
#   X3 <- rnorm(nsamp)
#   X4 <- rnorm(nsamp)
#   Z1 <- X1 + rnorm(nsamp,mean = 0,sd=1/4)
#   Z2 <- X2 + rnorm(nsamp,mean = 0,sd=1/4)
#   Z3 <- X3 + rnorm(nsamp,mean = 0,sd=1/4)
#   Z4 <- X4 + rnorm(nsamp,mean = 0,sd=1/4)
#
#   #put together into design matrix
#   X <- cbind(X1,X2,X3,X4,Z1,Z2,Z3,Z4)
#
#   #calculate Y
#   Y <- X%*%cbind(beta_true)+rnorm(nsamp)
#
#   #calculate variables shapley values
#   Yshap <- shapley(Y,X)
#   AICvectors[[i]]<- Yshap$AIC_shapleys
#   rsqvectors[[i]]<- Yshap$rsq_shapleys
# }
end_time = Sys.time()
sink(file="log.txt",append=TRUE,type="output")
cat("success!\n")
end_time
sink()

#post-simulation organizing#####################################################
#store shapley values into a data frame
AICshapleys <- finalMatrix[,1:8]
AICshapleys <- as.data.frame(AICshapleys)
AICshapleys$scorefn <- rep("AIC",nsim)
rsqshapleys <- finalMatrix[,9:16]
rsqshapleys <- as.data.frame(rsqshapleys)
rsqshapleys$scorefn <- rep("R_Squared",nsim)
shapleys <- dplyr::full_join(AICshapleys,rsqshapleys,by = c("X1", "X2", "X3",
                                                            "X4", "Z1", "Z2",
                                                            "Z3", "Z4",
                                                            "scorefn"))

#Store Results##################################################################
save(
  shapleys,
  file=output_file
)

#Log results####################################################################
sink(file="log.txt",append=TRUE,type="output")
cat("Simulation complete! Time elapsed:\n")
cat(end_time-start_time,"\n")
cat("Workspace saved at: ","\"",output_file,"\"\n")
sink()
