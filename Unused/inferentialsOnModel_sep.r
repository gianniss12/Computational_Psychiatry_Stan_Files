library(rstan)
library(parallel)
library(matlab)

setwd('/home/jonny/Documents/Internship-Gillan/Modeling-eLife')
load('stanfits_jonny_output_s1_Feb2_2020_v2.RData')
m = extract(fit, permuted=TRUE)
bySubj = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha1=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))

for( i in 1:length(unique(data$subj))){
    bySubj[i,2:6] = c(mean(m$betac[,i]),mean(m$alpha1[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]))
}

bySubj$ss = as.numeric(as.character(substr(bySubj$subj,2,5)))

sreps = T1_data[c("ss", "ocd", "gad", "T1_Comp", "T1_Gen", "T1_Obs", "T1_OCI_total", "T1_MCQ_total", "T1_DASS_total", "T1_sheehan_123", "T1_age")]
sreps$ss = as.numeric(as.character(sreps$ss))
comb_1 = merge(sreps, bySubj, by="ss")

#anxiety predicts change as a result of gas, p=.042
data=comb_2

summary(lm(beta1t~T1_Comp + T1_age, data=data))
summary(lm(alpha1~T1_Comp + T1_age, data=data))
summary(lm(betac~T1_Comp + T1_age, data=data))


summary(lm(beta1t~T1_Gen + T1_age, data=data))
summary(lm(alpha1~T1_Gen + T1_age, data=data))
summary(lm(betac~T1_Gen + T1_age, data=data))


summary(lm(beta1t~T1_Obs + T1_age, data=data))
summary(lm(alpha1~T1_Obs + T1_age, data=data))
summary(lm(betac~T1_Obs + T1_age, data=data))


summary(lm(beta1m~T1_Comp + T1_age, data=data))
summary(lm(beta1m~T1_Gen + T1_age, data=data))
summary(lm(beta1m~T1_Obs + T1_age, data=data))
summary(lm(beta1m~ocd + T1_age, data=data))
summary(lm(beta1m~gad + T1_age, data=data))
summary(lm(beta1m~ocd +T1_Comp+ T1_age, data=data))



load('stanfits_dps_sep_decay_T2.RData')
m = extract(fit, permuted=TRUE)
bySubj = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha1=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))

for( i in 1:length(unique(data$subj))){
    bySubj[i,2:6] = c(mean(m$betac[,i]),mean(m$alpha1[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]))
}

bySubj$ss = as.numeric(as.character(substr(bySubj$subj,2,5)))

sreps = T2_data[c("ss", "ocd", "gad", "T2_Comp", "T2_Gen", "T2_Obs", "T2_OCI_total", "T2_MCQ_total", "T2_DASS_total", "T2_sheehan_123", "T2_age")]
sreps$ss = as.numeric(as.character(sreps$ss))
comb_2 = merge(sreps, bySubj, by="ss")
