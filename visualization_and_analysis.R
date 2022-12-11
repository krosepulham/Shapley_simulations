library(shapleyAIC)
library(tidyverse)

#Trivial model##################################################################
load("Results/trivial.RData")

#Histograms
shapleys%>%
  pivot_longer(1:8,names_to = "variable",values_to = "shapley_value")%>%
  mutate(var_useful = grepl("X",variable))%>%
  ggplot(aes(x=shapley_value,fill=var_useful))+
  geom_histogram(position = "dodge")+
  facet_wrap(~scorefn,scales="free")+
  labs(
    title="Shapley Values for all slopes = 1",
    x = "Shapley Value",
    y = "Frequency",
    fill="Non-Zero slope?"
    )
ggsave(
  file="trivial-histogram.png",
  device="png",
  path="Figures",
  height = 4,
  width=8
)

#data frame for false positives
nsim <- nrow(filter(shapleys,scorefn=="AIC"))
AICshapleys <- filter(shapleys,scorefn=="AIC")
rsqshapleys <- filter(shapleys,scorefn=="R_Squared")
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

#summarize these two false-positive rates
apply(false_positives,MARGIN = 2,FUN=summary)
apply(false_positives,MARGIN=2,FUN=table)

#plot false positives for both schemes
ggplot(false_positives,aes(x=r_squared,y=AIC))+
  geom_jitter(height=0.05,width=0.05,alpha=0.2)+
  geom_abline(slope=1,intercept=0,linetype=2,color="red")+
  labs(
    title="False Positives for two methods",
    subtitle = "Each point is one simulation. Number of false positives on x and y axes for the two methods."
    )

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

#summarize these two separation ratios
apply(separation,MARGIN = 2,FUN=summary)

#plot these separation ratios, along with a y=x reference line to see when
#AIC shapley values achieve "Better separation"
ggplot(separation,aes(x=r_squared,y=AIC))+
  geom_point(alpha=0.3)+
  geom_abline(slope=1,intercept=0,linetype=2,color="red")+
  labs(
    title="\"Separation ratios\" for the two methods",
    subtitle="(Larger is better)"
  )


#one-half-trivial model#########################################################

load("Results/half-trivial.RData")

#Histograms
shapleys%>%
  pivot_longer(1:8,names_to = "variable",values_to = "shapley_value")%>%
  mutate(var_useful = grepl("X",variable))%>%
  ggplot(aes(x=shapley_value,fill=var_useful))+
  geom_histogram(position = "dodge")+
  facet_wrap(~scorefn,scales="free")+
  labs(
    title="Shapley Values for all slopes = 1/2",
    x = "Shapley Value",
    y = "Frequency",
    fill="Non-Zero slope?"
  )

#data frame for false positives
nsim <- nrow(filter(shapleys,scorefn=="AIC"))
AICshapleys <- filter(shapleys,scorefn=="AIC")
rsqshapleys <- filter(shapleys,scorefn=="R_Squared")
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

#summarize these two false-positive rates
apply(false_positives,MARGIN = 2,FUN=summary)
apply(false_positives,MARGIN=2,FUN=table)

#plot false positives for both schemes
ggplot(false_positives,aes(x=r_squared,y=AIC))+
  geom_jitter(height=0.05,width=0.05,alpha=0.2)+
  geom_abline(slope=1,intercept=0,linetype=2,color="red")+
  labs(
    title="False Positives for two methods",
    subtitle = "Each point is one simulation. Number of false positives on x and y axes for the two methods."
  )

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

#summarize these two separation ratios
apply(separation,MARGIN = 2,FUN=summary)

#plot these separation ratios, along with a y=x reference line to see when
#AIC shapley values achieve "Better separation"
ggplot(separation,aes(x=r_squared,y=AIC))+
  geom_point(alpha=0.3)+
  geom_abline(slope=1,intercept=0,linetype=2,color="red")+
  labs(
    title="\"Separation ratios\" for the two methods",
    subtitle="(Larger is better)"
  )

#1-1/2-1/3-1/4 model#########################################################

load("Results/one-half-third-fourth.RData")

#Histograms
shapleys%>%
  pivot_longer(1:8,names_to = "variable",values_to = "shapley_value")%>%
  mutate(var_useful = grepl("X",variable))%>%
  ggplot(aes(x=shapley_value,fill=var_useful))+
  geom_histogram(position = "dodge")+
  facet_wrap(~scorefn,scales="free")+
  labs(
    title="Shapley Values for all slopes = 1/2",
    x = "Shapley Value",
    y = "Frequency",
    fill="Non-Zero slope?"
  )

#data frame for false positives
nsim <- nrow(filter(shapleys,scorefn=="AIC"))
AICshapleys <- filter(shapleys,scorefn=="AIC")
rsqshapleys <- filter(shapleys,scorefn=="R_Squared")
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

#summarize these two false-positive rates
apply(false_positives,MARGIN = 2,FUN=summary)
apply(false_positives,MARGIN=2,FUN=table)

#plot false positives for both schemes
ggplot(false_positives,aes(x=r_squared,y=AIC))+
  geom_jitter(height=0.05,width=0.05,alpha=0.2)+
  geom_abline(slope=1,intercept=0,linetype=2,color="red")+
  labs(
    title="False Positives for two methods",
    subtitle = "Each point is one simulation. Number of false positives on x and y axes for the two methods."
  )

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

#summarize these two separation ratios
apply(separation,MARGIN = 2,FUN=summary)

#plot these separation ratios, along with a y=x reference line to see when
#AIC shapley values achieve "Better separation"
ggplot(separation,aes(x=r_squared,y=AIC))+
  geom_point(alpha=0.3)+
  geom_abline(slope=1,intercept=0,linetype=2,color="red")+
  labs(
    title="\"Separation ratios\" for the two methods",
    subtitle="(Larger is better)"
  )

