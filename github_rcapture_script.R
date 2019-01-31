rm(list=ls())
library(Rcapture)
cap_mat<-read.csv('github_rcapture_sample.csv') # subset of Jillian A. Hudgins' collated turtle photo ID repository for the Republic of Maldives, each row is a different turtle and each column is a pre-specified minor capture interval (1 month in this case)
rownames(cap_mat)<-NULL #RCapture requires NULL row names

# Separate into 6-month major capture intervals
# dfreq specifies whether the last column of the matrix includes total counts, which is not the case for our input matrix
# drop (default) removes turtles without any sightings
# vt specifies the capture interval structure (8 blocks of 6 months in our case) and must add up to the total size of the matrix
reduced<-periodhist(cap_mat, dfreq=FALSE,vt=c(rep(6,8)),drop=TRUE)
hist(reduced[,ncol(reduced)], xlab="number of sightings", main=NULL) #plot raw frequencies

#openp calculates the best fitting open population model, (notice that reduced has a frequency column so dfreq is now TRUE)
#m specifies model type, 'ep' or equal probability of capture vs. heterogeneous probability ('up') across major intervals 
#neg defaults to true and keeps survival on (0,1) and births positive
open_model_ep<-openp(reduced, dfreq=TRUE, m="ep", neg=TRUE)
open_model_up<-openp(reduced, dfreq=TRUE, m="up", neg=TRUE)

open_model_ep$trap.param # trap effect fitted parameter
open_model_ep$model.fit # model without trap effect
open_model_ep$trap.fit # model with trap effect

open_model_up$trap.param # trap effect fitted parameter
open_model_up$model.fit # model without trap effect
open_model_up$trap.fit # model with trap effect

# best open model based on AIC us unequal probability model without a trap effect (AIC=149.27)

closed_model<-closedp(reduced, dfreq=TRUE, neg=TRUE)
closed_model$results # AIC is higher for all closed models, so best model is open_model_up


## Best model results ##
open_model_up$survivals # survival rate
open_model_up$Ntot # total pop size ever recorded
open_model_up$N # pop size per major capture interval
open_model_up$birth # births
open_model_up$capture.prob # probability of capture
