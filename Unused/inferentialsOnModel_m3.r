library(rstan)
library(parallel)
library(matlab)

setwd('/Users/clairegillan/OneDrive - TCDUD.onmicrosoft.com/Gillan Lab Resources/Projects/Model_Test_Jonny/')
load('S2_results_withdecay_sepBetas_noLambda_nodecay_noRescaling.RData')
m = extract(fit, permuted=TRUE)
bySubj = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))),alpha1=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))

for( i in 1:length(unique(data$subj))){
    bySubj[i,2:6] = c(mean(m$betac[,i]),mean(m$alpha1[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]))
}

bySubj$ss = as.character(bySubj$subj)

# set up fot s2 data:
sr_data = read.csv('twostepForModeling_s2.csv')

sr_data$ss = sr_data$subject
sr_data = sr_data[!duplicated(sr_data$ss),]
comb_1 = merge(sr_data, bySubj, by="ss")


data=comb_1
summary(lm(beta1m~scale(ocd) + scale(age) + scale(iq) + gender, data=data))
summary(lm(beta1m~scale(alcoholAddiction) + scale(age) + scale(iq) + gender, data=data))
summary(lm(beta1m~scale(eatingDisorders) + scale(age) + scale(iq) + gender, data=data))
summary(lm(beta1m~scale(impulsivity) + scale(age) + scale(iq) + gender, data=data))

# linking to the transdiagnostic dimension for s2 data:

load('efa_results_scaled_f3.RData')
comb_1$ss = as.factor(paste("s",comb_1$ss,sep=""))
comb_2 = merge(comb_1, small_s2, by="ss")

summary(lm(beta1m~scale(Factor2) + agez + iqz + gender.y, data=comb_2))
