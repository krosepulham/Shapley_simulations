#setup##########################################################################
rm(list = ls())
set.seed(232826833)
AICvectors <- list()
rsqvectors <- list()
nsim  <- 1000
nsamp <- 50
beta_true <- c(1,1,1,1,0,0,0,0)
start_time = Sys.time()
output_file <- "Results/x-correlated.RData"

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
for(i in 1:nsim){
  #randomly generate explanatory covariates
  meanvec <- c(0,0,0,0)
  rho <- 0.7
  varcovar <- rbind(
    c(1  ,rho,rho,-rho),
    c(rho,1  ,rho,-rho),
    c(rho,rho,1  ,-rho),
    c(-rho,-rho,-rho,1  )
  )
  X <- MASS::mvrnorm(n=50,mu=meanvec,Sigma=varcovar)

  X1 <- X[,1]
  X2 <- X[,2]
  X3 <- X[,3]
  X4 <- X[,4]
  Z1 <- rnorm(nsamp)
  Z2 <- rnorm(nsamp)
  Z3 <- rnorm(nsamp)
  Z4 <- rnorm(nsamp)

  #put together into design matrix
  X <- cbind(X1,X2,X3,X4,Z1,Z2,Z3,Z4)

  #calculate Y
  Y <- X%*%cbind(beta_true)+rnorm(nsamp)

  #calculate variables shapley values
  Yshap <- shapley(Y,X)
  AICvectors[[i]]<- Yshap$AIC_shapleys
  rsqvectors[[i]]<- Yshap$rsq_shapleys
}
end_time = Sys.time()
sink(file="log.txt",append=TRUE,type="output")
cat("success!\n")
end_time
sink()

#post-simulation organizing#####################################################
#store shapley values into a data frame
AICshapleys <- matrix(unlist(AICvectors),nrow=nsim,byrow=TRUE)
colnames(AICshapleys) <- names(AICvectors[[1]])
AICshapleys <- as.data.frame(AICshapleys)
AICshapleys$scorefn <- rep("AIC",nsim)
rsqshapleys <- matrix(unlist(rsqvectors),nrow=nsim,byrow=TRUE)
colnames(rsqshapleys) <- names(rsqvectors[[1]])
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
