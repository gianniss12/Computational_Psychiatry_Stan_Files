library(rstan)
library(parallel)
library(matlab)
library(shinystan)

## Load in Stan Fitted Model ##################

setwd('C:/Users/jgior/Desktop/Internship/Modeling-eLife/')
load('stan_fits/s1_thesis/JG_s1_sepAlphas_sepStages_qt2-qt1_zeroLambda_noDecay_March8_2020_mistake.RData')


### Check R-Hat Values ######################
fit_summary = summary(fit)
print(fit_summary$summary)

#Shinystan
my_shinystan = as.shinystan(fit)
launch_shinystan(my_shinystan)


### Derive subject means for Betas, Alpha #######################
m = extract(fit, permuted=TRUE)

### For Single Alpha ##########################
bySubj = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha1=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
 
for( i in 1:length(unique(data$subj))){
     bySubj[i,2:6] = c(mean(m$betac[,i]),mean(m$alpha1[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]))
}

### For Multiple Alphas, Both Stages ##########################

bySubj = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha_p=rep(0,(length(unique(data$subj)))), alpha_n=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
for( i in 1:length(unique(data$subj))){
bySubj[i,2:7] = c(mean(m$betac[,i]),mean(m$alpha_p[,i]), mean(m$alpha_n[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]))
}
summary(bySubj)

## Checking 2.5% Quantile
bySubj_2.5 = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha_p=rep(0,(length(unique(data$subj)))), alpha_n=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
for( i in 1:length(unique(data$subj))){
  bySubj_2.5[i,2:7] = c(quantile(m$betac[,i], probs = c(.025)),quantile(m$alpha_p[,i], probs = c(.025)), quantile(m$alpha_n[,i], probs = c(.025)), quantile(m$beta1m[,i], probs = c(.025)),quantile(m$beta1t[,i], probs = c(.025)), quantile(m$beta2[,i], probs = c(.025)))
}
summary(bySubj_2.5)

##Checking 97.5% Quantile
bySubj_97.5 = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha_p=rep(0,(length(unique(data$subj)))), alpha1_n=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
for( i in 1:length(unique(data$subj))){
  bySubj_97.5[i,2:7] = c(quantile(m$betac[,i], probs = c(.975)),quantile(m$alpha_p[,i], probs = c(.975)), quantile(m$alpha_n[,i], probs = c(.975)), quantile(m$beta1m[,i], probs = c(.975)),quantile(m$beta1t[,i], probs = c(.975)), quantile(m$beta2[,i], probs = c(.975)))
}
summary(bySubj_97.5)

### For Multiple Alphas, Seperate Stages ############################

## Calculate Mean
bySubj = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha1_p=rep(0,(length(unique(data$subj)))), alpha1_n=rep(0,(length(unique(data$subj)))), alpha2_p=rep(0,(length(unique(data$subj)))), alpha2_n=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
for( i in 1:length(unique(data$subj))){
  bySubj[i,2:9] = c(mean(m$betac[,i]),mean(m$alpha1_p[,i]), mean(m$alpha1_n[,i]), mean(m$alpha2_p[,i]), mean(m$alpha2_n[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]))
}
summary(bySubj)

## Checking 2.5% Quantile
bySubj_2.5 = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha1_p=rep(0,(length(unique(data$subj)))), alpha1_n=rep(0,(length(unique(data$subj)))), alpha2_p=rep(0,(length(unique(data$subj)))), alpha2_n=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
for( i in 1:length(unique(data$subj))){
  bySubj_2.5[i,2:9] = c(quantile(m$betac[,i], probs = c(.025)),quantile(m$alpha1_p[,i], probs = c(.025)), quantile(m$alpha1_n[,i], probs = c(.025)), quantile(m$alpha2_p[,i], probs = c(.025)), quantile(m$alpha2_n[,i], probs = c(.025)), quantile(m$beta1m[,i], probs = c(.025)),quantile(m$beta1t[,i], probs = c(.025)), quantile(m$beta2[,i], probs = c(.025)))
}
summary(bySubj_2.5)


##Checking 97.5% Quantile
bySubj_97.5 = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha1_p=rep(0,(length(unique(data$subj)))), alpha1_n=rep(0,(length(unique(data$subj)))), alpha2_p=rep(0,(length(unique(data$subj)))), alpha2_n=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
for( i in 1:length(unique(data$subj))){
  bySubj_97.5[i,2:9] = c(quantile(m$betac[,i], probs = c(.975)),quantile(m$alpha1_p[,i], probs = c(.975)), quantile(m$alpha1_n[,i], probs = c(.975)), quantile(m$alpha2_p[,i], probs = c(.975)), quantile(m$alpha2_n[,i], probs = c(.975)), quantile(m$beta1m[,i], probs = c(.975)),quantile(m$beta1t[,i], probs = c(.975)), quantile(m$beta2[,i], probs = c(.975)))
}
summary(bySubj_97.5)


## Load in Psychiatric Scores for S1 data #########################

scores = read.csv('data/self_report_study1.csv', header=TRUE)
sreps = scores[c("subj", "age", "iq", "gender", "sds_total", "stai_total", "oci_total")]

## Merge Psychiatric scores with Subject fittings
comb = merge(sreps, bySubj, by="subj")

###Analysis - Model-Based Beta

summary(lm(beta1m ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm(beta1m ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(beta1m ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))


### Analysis for 2 Learning Rates ###
## Analysis - Alpha_p

summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))

summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))

## Analysis - Alpha_n

summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))

## Analysis - [alpha_n - alpha_p]

summary(lm((alpha_n-alpha_p) ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm((alpha_n-alpha_p) ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm((alpha_n-alpha_p) ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))


### Analysis for 4 Learning Rates (2 LRs per Step) ###
## Analysis - Alpha1_p

summary(lm(alpha1_p ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm(alpha1_p ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(alpha1_p ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))

## Analysis - Alpha1_n

summary(lm(alpha1_n ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm(alpha1_n ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(alpha1_n ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))

## Analysis - [alpha1_n - alpha1_n]

summary(lm((alpha1_n-alpha1_p) ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm((alpha1_n-alpha1_p) ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm((alpha1_n-alpha1_p) ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))

## Analysis - Alpha2_p

summary(lm(alpha2_p ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm(alpha2_p ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(alpha2_p ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))

## Analysis - Alpha2_n

summary(lm(alpha2_n ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm(alpha2_n ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(alpha2_n ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))

## Analysis - [alpha2_n - alpha2_p]

summary(lm((alpha2_n-alpha2_p) ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm((alpha2_n-alpha2_p) ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm((alpha2_n-alpha2_p) ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))


