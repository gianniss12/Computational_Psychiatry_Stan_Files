library(rstan)
library(parallel)
library(matlab)

setwd('/Users/clairegillan/Documents/Model_Test_Jonny/')
load('TScode_elife_s1results.RData')
m = extract(fit, permuted=TRUE)
bySubj = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))),
lambda=rep(0,(length(unique(data$subj)))),alpha1=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))

for( i in 1:length(unique(data$subj))){
    bySubj[i,2:7] = c(mean(m$betac[,i]),mean(m$lambda[,i]),mean(m$alpha1[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]))
}

bySubj$ss = as.character(bySubj$subj)

sr_data = read.csv('twostepForModeling_s1.csv')

sr_data$ss = sr_data$subject
sr_data = sr_data[!duplicated(sr_data$ss),]
comb_1 = merge(sr_data, bySubj, by="ss")


data=comb_1
summary(lm(beta1m~scale(ocd) + scale(age) + scale(iq) + gender, data=data))
summary(lm(beta1m~scale(traitAnxiety) + scale(age) + scale(iq) + gender, data=data))
summary(lm(beta1m~scale(depression) + scale(age) + scale(iq) + gender, data=data))

summary(lm(alpha1~scale(depression) + scale(age) + scale(iq) + gender, data=data))



