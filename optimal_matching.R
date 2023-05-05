#Here we consider matching participants from a group where group is defined by
#type of exercise performed by the participants to, a control group
#the control group is defined by all the participants who enrolled in the study
#
#The control group is chosen as such because we have data for all the participants
#before they started doing any of the recommended exercises.
#
#We can match each treatment group to all the participants.
#The control group and the treatment group will be matched based on:
#RHR, HRV, Sleep Score, Hours of sleep, sleep efficacy, hours of sleep and
#respiration rate
#Outcomes before and after the exercises will be chosen for determining 
#treatment effect
#-------------------------------------------------------------------------------
#Preparing dataset
library(readxl)
library(dplyr)
library(tidyverse)
filename <- read_excel("BWPilot_CombinedData_20210803_fordryad_addvars_cleaned_noround2_v3.xlsx")
filename[filename=='.'] <- NA
data_set <- filename[, c(1,2,15,16,17,18,19,20)]


data_set <- data_set %>% drop_na()

data_set$`Days from  Round1 Day1` <- as.numeric(data_set$`Days from  Round1 Day1`)

data_set$RHR <- as.numeric(data_set$RHR)
data_set$HRV <- as.numeric(data_set$HRV)
data_set$`Sleep Score` <- as.numeric(data_set$`Sleep Score`)
data_set$`Hours of Sleep` <- as.numeric(data_set$`Hours of Sleep`)
data_set$`Sleep Efficiency` <- as.numeric(data_set$`Sleep Efficiency`)
data_set$`Respiration Rate` <- as.numeric(data_set$`Respiration Rate`)

t.first <- filename[match(unique(data_set$ScrSubjectID), data_set$ScrSubjectID),c(3)]
t.sec <- ifelse(t.first$`Round 1 Exercise` == "Box Breathing", 1,
                ifelse(t.first$`Round 1 Exercise` == "Slow Breathing", 2,
                       ifelse(t.first$`Round 1 Exercise` == "Super Oxygenation", 3, 4)) 
                )

a <- rep(0, 113)

treatment <- c(rbind(a, t.sec))

treatment <- treatment[-122]
treatment <- treatment[-137]

data_set <- data_set %>% 
  mutate(pos = `Days from  Round1 Day1`>0)

data_set_matching <- data_set %>% group_by(ScrSubjectID, pos) %>%
  summarise_all(mean)

data_set_matching <- data.frame(data_set_matching, treatment)
#-------------------------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(MASS)

subset <- data_set_matching[data_set_matching$treatment %in% c(0,1), ]  
datatemp=subset

subset$Hours.of.Sleep
subset$Sleep.Efficiency
subset$Respiration.Rate

#Score model
propscore.model=glm(pos ~ RHR + HRV + Sleep.Score +
                    Sleep.Efficiency + Respiration.Rate,
                    family=binomial,x=TRUE,y=TRUE,
                    data=subset)
# #-----------------------------------------------------------------------------
datatemp$treated=propscore.model$y
datatemp$treatment=datatemp$treated

# Use the caret package to include all categories of categorical variables (i.e.,
# do not leave out one category) in X matrix
library(caret)
dmy=dummyVars(propscore.model$formula,data=datatemp)
Xmat=data.frame(predict(dmy,newdata=datatemp))
# Matrix of covariates to include in the Mahalanobis distance, based on prognostically
#important variables
Xmatmahal=Xmat
treated=datatemp$treated
datatemp$logit.ps=predict(propscore.model)
# Use Hansen (2009)’s rule for removing subjects who lack overlap
# from Notes 6
logit.propscore=predict(propscore.model)
pooled.sd.logit.propscore=sqrt(var(logit.propscore[datatemp$treatment==1])/2+var(logit.propscore[datatemp$treatment==0])/2)
min.treated.logit.propscore=min(logit.propscore[datatemp$treatment==1])
max.control.logit.propscore=max(logit.propscore[datatemp$treatment==0])
# How many treated and control subjects lack overlap by Hansen's criterion
no.treated.lack.overlap=sum(logit.propscore[datatemp$treatment==1]>(max.control.logit.propscore+.5*pooled.sd.logit.propscore))
no.control.lack.overlap=sum(logit.propscore[datatemp$treatment==0]<(min.treated.logit.propscore-.5*pooled.sd.logit.propscore))
# If there are subjects who lack overlap, remove them from the datatemp dataset
datatemp.original=datatemp
datatemp.full=datatemp
Xmat.original=Xmat
Xmat.full=Xmat
if(no.treated.lack.overlap+no.control.lack.overlap>0){
  which.remove=which((logit.propscore>(max.control.logit.propscore+.5*pooled.sd.logit.propscore))|(logit.propscore<(min.treated.logit.propscore-.5*pooled.sd.logit.propscore)))
  datatemp=datatemp[-which.remove,]
  datatemp.full=rbind(datatemp,datatemp.original[which.remove,])
  Xmat=Xmat[-which.remove,]
  Xmat.full=rbind(Xmat,Xmat.original[which.remove,])
  Xmatmahal=Xmatmahal[-which.remove,]
}
rownames(datatemp)=seq(1,nrow(datatemp),1)

smahal=
  function(z,X){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z)
    for (j in 1:k) X[,j]<-rank(X[,j])
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    out<-matrix(NA,m,n-m)
    Xc<-X[z==0,]
    Xt<-X[z==1,]
    9
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
    icov<-ginv(cv)
    for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    out
  }

addcaliper=function(dmat,z,logitp,calipersd=.5,penalty=1000){
  # Pooled within group standard devation
  sd.logitp=sqrt((sd(logitp[z==1])^2+sd(logitp[z==0])^2)/2)
  adif=abs(outer(logitp[z==1],logitp[z==0],"-"))
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat=dmat+adif*penalty
  dmat
}

# Rank based Mahalanobis distance
distmat=smahal(datatemp$treated,Xmatmahal)
# Add caliper
distmat=addcaliper(distmat,datatemp$treated,datatemp$logit.ps,calipersd=.5)
# Label the rows and columns of the distance matrix by the rownames in datatemp
rownames(distmat)=rownames(datatemp)[datatemp$treated==1]
colnames(distmat)=rownames(datatemp)[datatemp$treated==0]
#-----------------#-----------------#-----------------#-----------------
library(optmatch)
matchvec=fullmatch(distmat,data=datatemp) 
datatemp$matchvec=matchvec 
# Number the strata
matchedset.index=substr(matchvec,start=3,stop=10)
matchedset.index.numeric=as.numeric(matchedset.index)
# Have matchedset.index.numeric.full append 0 to matchedset index.numeric for
# the removed subjects
if(no.control.lack.overlap+no.treated.lack.overlap==0){
  matchedset.index.numeric.full=matchedset.index.numeric
}
if(no.control.lack.overlap+no.treated.lack.overlap>0){
  matchedset.index.numeric.full=c(matchedset.index.numeric,rep(0,no.control.lack.overlap+no.treated.lack.overlap))
}
# Calculate standardized difference before and after a full match
# Calculate standardized difference before and after a full match
# Drop observations with missing values from the calculations
# stratum.myindex should contain strata for each subject, 0 means a unit was not
# matched
# Use harmonic mean weights
standardized.diff.harmonic.func=function(x,treatment,stratum.myindex,missing=rep(0,length(x))){
  xtreated=x[treatment==1 & missing==0];
  xcontrol=x[treatment==0 & missing==0];
  var.xtreated=var(xtreated);
  var.xcontrol=var(xcontrol);
  combinedsd=sqrt(.5*(var.xtreated+var.xcontrol));
  std.diff.before.matching=(mean(xtreated)-mean(xcontrol))/combinedsd;
  nostratum=length(unique(stratum.myindex))-1*max(stratum.myindex==0);
  if(max(stratum.myindex==0)==0){
    stratumlist=sort(unique(stratum.myindex))
  }
  if(max(stratum.myindex==0)==1){
    templist=sort(unique(stratum.myindex))
    stratumlist=templist[-1]
  }
  diff.in.stratum=rep(0,nostratum);
  number.in.stratum=rep(0,nostratum);
  harmonic.weight=rep(0,nostratum)
  for(i in 1:nostratum){
    if(sum(stratum.myindex==stratumlist[i] & treatment==1 & missing==0)==0 | sum(stratum.myindex==stratumlist[i] & treatment==0 & missing==0)==0){
      number.in.stratum[i]=0
    }
    if(sum(stratum.myindex==stratumlist[i] & treatment==1 & missing==0)>0 & sum(stratum.myindex==stratumlist[i] & treatment==0 & missing==0)>0){
      diff.in.stratum[i]=mean(x[stratum.myindex==stratumlist[i] & treatment==1 & missing==0])-mean(x[stratum.myindex==stratumlist[i] & treatment==0 & missing==0]);
      number.in.stratum[i]=sum(stratum.myindex==stratumlist[i])
      harmonic.weight[i]=1/(.5/sum(stratum.myindex==stratumlist[i] & treatment==1)+.5/sum(stratum.myindex==stratumlist[i] & treatment==0))
    }
  }
  std.diff.after.matching=(sum(harmonic.weight*diff.in.stratum)/sum(harmonic.weight))/combinedsd;
  list(std.diff.before.matching=std.diff.before.matching,std.diff.after.matching=std.diff.after.matching);
}

# Also compute balance on logit propensity score
# Xmatmahal$logit.ps=datatemp$logit.ps
# Calculate the standardized differences
std.diff.before=rep(0,ncol(Xmatmahal));
std.diff.after=rep(0,ncol(Xmatmahal));
names(std.diff.before)=names(Xmatmahal[1,]);
names(std.diff.after)=names(Xmatmahal[1,]);
for(i in 1:ncol(Xmatmahal)){
  missing.temp=is.na(Xmatmahal[,i])
  temp.stand.diff=standardized.diff.harmonic.func(Xmatmahal[,i],datatemp$treated,matchedset.index.numeric,missing.temp);
  std.diff.before[i]=temp.stand.diff$std.diff.before.matching;
  std.diff.after[i]=temp.stand.diff$std.diff.after.matching;
}
# Rename std.diff.before and std.diff.after to shorter names sd.bf and sd.af
# and use digits option to be able to columns of std.diff.before and
# std.diff.after in one row
sd.bf=std.diff.before
sd.af=std.diff.after
options(digits=2)
cbind(sd.bf,sd.af)

covariates=names(sd.bf)

plot.dataframe2=data.frame(stand.diff=c(sd.bf,sd.af),covariates=rep(covariates,2),type=c(rep("Before",length(covariates)),rep("After",length(covariates))))
ggplot(plot.dataframe2,aes(x=stand.diff,y=covariates))+geom_point(size=5,aes(shape=factor(type)))+scale_shape_manual(values=c(4,1))+geom_vline(xintercept=c(-0.2,0.2),lty=2)+
  theme(axis.text.y = element_text(size=6)) 

ggsave(
  filename = "l5.png",
  device = "png", width = 6, height = 4)

# #Point estimate for optimal matching and Confidence interval.
# reg.formula=update(propscore.model$formula,outcome~treated+matchvec+.)
# matched.reg.model=lm(reg.formula,data=datatemp)
# summary(matched.reg.model)
# # Point estimate of treatment effect
# coef(matched.reg.model)[2]
# # Confidence interval
# confint(matched.reg.model)[2,]


