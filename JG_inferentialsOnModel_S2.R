library(rstan)
library(parallel)
library(matlab)

## Load in Stan Fitted Model

setwd('C:/Users/jgior/Desktop/Internship/Modeling-eLife/')
load('stan_fits/stanfits_jonny_output_s2_Feb8_1500iter_2020.RData')

## Check R-Hat Values

fit_summary = summary(fit)
print(fit_summary$summary)

## Derive subject means for Betas, Alpha
m = extract(fit, permuted=TRUE)

# For Single Alpha
# bySubj = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha1=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
# for( i in 1:length(unique(data$subj))){
#     bySubj[i,2:6] = c(mean(m$betac[,i]),mean(m$alpha1[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]))
# }

# For Multiple Alphas
bySubj = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha_p=rep(0,(length(unique(data$subj)))), alpha_n=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
for( i in 1:length(unique(data$subj))){
  bySubj[i,2:7] = c(mean(m$betac[,i]),mean(m$alpha_p[,i]), mean(m$alpha_n[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]))
}


## Load in Psychiatric Scores

scores = read.csv('data/self_report_study2.csv', header=TRUE)
sreps = scores[c("subj", "age","iq", "gender", "sds_total","stai_total","oci_total", 
                 "lsas_total", "bis_total", "scz_total", "aes_total", "eat_total",
                 "audit_total", "Factor1","Factor2","Factor3")]

## Merge Psychiatric scores with Subject fittings
comb = merge(sreps, bySubj, by="subj")

## Analysis

## GLM for Factors

summary(lm(beta1m ~ scale(iq) + scale(age) + gender + scale(Factor1), data=comb))
summary(lm(beta1m ~ scale(iq) + scale(age) + gender + scale(Factor2), data=comb))
summary(lm(beta1m ~ scale(iq) + scale(age) + gender + scale(Factor3), data=comb))


## GLM for Scales

summary(lm(beta1m ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm(beta1m ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(beta1m ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))