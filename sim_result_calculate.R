#sim_result_calculate() function ###############################################
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
#function that takes a result file name and does the following:
#>plots side-by-side histograms of combined shapley value pools and ggsave()
#>false positive rates, plot and ggsave(), create tables and write to .txt
#>correlation ratios, plot and ggsave(), create tables and write to .txt
sim_result_calculate<-function(filename){
  # Load data and calculate false positive and separation ratios ###############
  load(str_glue("Sim_results/",filename))
  setting <- str_split_1(filename,pattern = ".RData")[1]

  #data frame for false positives
  nsim <- nrow(dplyr::filter(shapleys,scorefn=="AIC"))
  AICshapleys <- dplyr::filter(shapleys,scorefn=="AIC")
  rsqshapleys <- dplyr::filter(shapleys,scorefn=="R_Squared")
  false_positives <- data.frame(
    AIC = rep(NA,nsim),
    r_squared = rep(NA,nsim)
  )

  #for loop for AIC false positives
  for(i in 1:nsim){
    smallest_X <- min(AICshapleys[i,1:4])
    false_positives$AIC[i]<-sum(AICshapleys[i,5:8]>smallest_X)
  }

  #for loop for r-squared false positives
  for(i in 1:nsim){
    smallest_X <- min(rsqshapleys[i,1:4])
    false_positives$r_squared[i]<-sum(rsqshapleys[i,5:8]>smallest_X)
  }

  #Separation Ratios
  separation <- data.frame(
    AIC = rep(NA,nsim),
    r_squared = rep(NA,nsim)
  )

  #for loop for AIC false positives
  for(i in 1:nsim){
    smallest_X <- min(AICshapleys[i,1:4])
    largest_Z  <- max(AICshapleys[i,5:8])
    smallest   <- min(AICshapleys[i,1:8])
    largest    <- max(AICshapleys[i,1:8])
    separation$AIC[i]<- (smallest_X-largest_Z)/(largest-smallest)
  }

  #for loop for r-squared false positives
  for(i in 1:nsim){
    smallest_X <- min(rsqshapleys[i,1:4])
    largest_Z  <- max(rsqshapleys[i,5:8])
    smallest   <- min(rsqshapleys[i,1:8])
    largest    <- max(rsqshapleys[i,1:8])
    separation$r_squared [i]<- (smallest_X-largest_Z)/(largest-smallest)
  }

  # create plots ###############################################################

  #histogram
  phist <- shapleys%>%
    pivot_longer(1:8,names_to = "variable",values_to = "shapley_value")%>%
    mutate(var_useful = grepl("X",variable))%>%
    ggplot(aes(x=shapley_value,fill=var_useful))+
    geom_histogram(position = "dodge",bins = 30)+
    facet_wrap(~scorefn,scales="free")+
    labs(
      title="Shapley Values for all slopes = 1",
      x = "Shapley Value",
      y = "Frequency",
      fill="Non-Zero slope?"
    )
  ggsave(
    filename = str_glue("Figures/",setting,"_histogram.png"),
    plot = phist,
    width = 6,
    height = 3,
    units="in"
  )

  #plot false positives for both schemes
  p1 <- ggplot(false_positives,aes(x=r_squared,y=AIC))+
    geom_jitter(height=0.05,width=0.05,alpha=0.2)+
    geom_abline(slope=1,intercept=0,linetype=2,color="red")+
    labs(
      title="False Positives for two methods",
      subtitle = "Each point is one simulation. Number of false positives on x and y axes for the two methods."
    )
  ggsave(
    filename = str_glue("Figures/",setting,"_false-positives.png"),
    plot = p1,
    width = 4,
    height = 4,
    units="in"
  )

  #AIC shapley values achieve "Better separation"
  p2 <- ggplot(separation,aes(x=r_squared,y=AIC))+
    geom_point(alpha=0.3)+
    geom_abline(slope=1,intercept=0,linetype=2,color="red")+
    labs(
      title="\"Separation ratios\" for the two methods",
      subtitle="(Larger is better)"
    )
  ggsave(
    filename = str_glue("Figures/",setting,"_separation-ratios.png"),
    plot = p2,
    width = 4,
    height = 4,
    units="in"
  )

  #Write summary stats and tables to a text file################################
  cat("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
  pstring <- str_glue("Results for ",setting,":\n")
  print(pstring)


  #summarize these two false-positive rates
  cat("False Positive rates:\n")
  print(apply(false_positives,MARGIN = 2,FUN=summary))
  print(t(apply(false_positives,MARGIN=2,FUN=table)))
  #summarize these two separation ratios
  cat("\n\nSeparation Ratios\n")
  print(apply(separation,MARGIN = 2,FUN=summary))


}

#use function on results
if (file.exists("Summaries_and_tables.txt")){
  file.remove("Summaries_and_tables.txt")
}
file.create("Summaries_and_tables.txt")
sink(file="Summaries_and_tables.txt",append = TRUE)
sim_result_calculate("trivial.RData")
sim_result_calculate("half-trivial.RData")
sim_result_calculate("one-half-third-fourth.RData")
sim_result_calculate("x-correlated.RData")
sim_result_calculate("x-z-correlated.RData")
sink()
