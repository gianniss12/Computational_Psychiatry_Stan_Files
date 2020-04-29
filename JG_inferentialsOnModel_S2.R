library(rstan)
library(parallel)
library(matlab)
library(shinystan)
library(ggplot2)

## Load in Stan Fitted Model ######################

setwd('C:/Users/jgior/Desktop/Internship/Modeling-eLife/')
load('stan_fits/s2_thesis/JG_s2_sepAlphas_r-qt1_oneLambda_noDecay_April15_2020.RData')



### Check R-Hat Values ######################
fit_summary = summary(fit)
print(fit_summary$summary)

#Shinystan
my_shinystan = as.shinystan(fit)
launch_shinystan(my_shinystan)



### Derive subject means for Betas, Alpha #######################
m = extract(fit, permuted=TRUE)

### For Single Alpha ##########################

#Checking Mean
bySubj = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha1=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
for( i in 1:length(unique(data$subj))){
  bySubj[i,2:6] = c(mean(m$betac[,i]),mean(m$alpha1[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]))
}
summary(bySubj)

## Checking 2.5% Quantile
bySubj_2.5 = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha1=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
for( i in 1:length(unique(data$subj))){
  bySubj_2.5[i,2:6] = c(quantile(m$betac[,i], probs = c(.025)),quantile(m$alpha1[,i], probs = c(.025)), quantile(m$beta1m[,i], probs = c(.025)),quantile(m$beta1t[,i], probs = c(.025)), quantile(m$beta2[,i], probs = c(.025)))
}
summary(bySubj_2.5)

##Checking 97.5% Quantile
bySubj_97.5 = data.frame(subj=unique(data$subj), betac=rep(0,(length(unique(data$subj)))), alpha1=rep(0,(length(unique(data$subj)))), beta1m=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))),beta2=rep(0,(length(unique(data$subj)))))
for( i in 1:length(unique(data$subj))){
  bySubj_97.5[i,2:6] = c(quantile(m$betac[,i], probs = c(.975)),quantile(m$alpha1[,i], probs = c(.975)), quantile(m$beta1m[,i], probs = c(.975)),quantile(m$beta1t[,i], probs = c(.975)), quantile(m$beta2[,i], probs = c(.975)))
}
summary(bySubj_97.5)

# Histogram
hist(bySubj$alpha1, xlab="Mean Learning Rate", ylab="Count", main="Histogram of Mean Learning Rates Across Subjects")





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


## Load in Psychiatric Scores #####################

scores = read.csv('data/self_report_study2.csv', header=TRUE)
sreps = scores[c("subj", "age","iq", "gender", "sds_total","stai_total","oci_total", 
                 "lsas_total", "bis_total", "scz_total", "aes_total", "eat_total",
                 "audit_total", "Factor1","Factor2","Factor3")]

## Merge Psychiatric scores with Subject fittings
comb = merge(sreps, bySubj, by="subj")

### Analysis ###########################


symptoms = c("sds","stai","oci", "lsas", "bis", "scz", "aes", "eat", "audit")
coeff = rep(0,9)

symp_p = data.frame(symptoms, coeff)
symp_n = data.frame(symptoms, coeff)


factors = c("Factor1", "Factor2", "Factor3")
coeff = rep(0,3)
fact = data.frame(factors, coeff)

### Analysis for 1 Learning Rate

summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(lsas_total), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(bis_total), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(scz_total), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(aes_total), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(eat_total), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(audit_total), data=comb))



summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(Factor1), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(Factor2), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(Factor3), data=comb))
summary(lm(alpha1 ~ scale(iq) + scale(age) + gender + scale(Factor1) + scale(Factor2) + scale(Factor3), data=comb))



### Analysis for 2 Learning Rates ###
## Analysis - Alpha_p

#Create list of scales, zerod coeff
symptoms = c("sds","stai","oci", "lsas", "bis", "scz", "aes", "eat", "audit")

coeff = rep(0,9)
coeff[1] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))$coefficients[5,1]
coeff[2] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))$coefficients[5,1]
coeff[3] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))$coefficients[5,1]
coeff[4] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(lsas_total), data=comb))$coefficients[5,1]
coeff[5] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(bis_total), data=comb))$coefficients[5,1]
coeff[6] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(scz_total), data=comb))$coefficients[5,1]
coeff[7] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(aes_total), data=comb))$coefficients[5,1]
coeff[8] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(eat_total), data=comb))$coefficients[5,1]
coeff[9] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(audit_total), data=comb))$coefficients[5,1]
coeff = coeff*100
symp_p = data.frame(symptoms, coeff)

coeff = rep(0,3)
coeff[1] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(Factor1), data=comb))$coefficients[5,1]
coeff[2] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(Factor2), data=comb))$coefficients[5,1]
coeff[3] = summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(Factor3), data=comb))$coefficients[5,1]
coeff = coeff*100
fact_p = data.frame(factors, coeff)


##########

summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(lsas_total), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(bis_total), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(scz_total), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(aes_total), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(eat_total), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(audit_total), data=comb))



summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(Factor1), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(Factor2), data=comb))
summary(lm(alpha_p ~ scale(iq) + scale(age) + gender + scale(Factor3), data=comb))



## Analysis - Alpha_n

coeff = rep(0,9)
coeff[1] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))$coefficients[5,1]
coeff[2] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))$coefficients[5,1]
coeff[3] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))$coefficients[5,1]
coeff[4] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(lsas_total), data=comb))$coefficients[5,1]
coeff[5] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(bis_total), data=comb))$coefficients[5,1]
coeff[6] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(scz_total), data=comb))$coefficients[5,1]
coeff[7] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(aes_total), data=comb))$coefficients[5,1]
coeff[8] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(eat_total), data=comb))$coefficients[5,1]
coeff[9] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(audit_total), data=comb))$coefficients[5,1]
coeff = coeff*100
symp_n = data.frame(symptoms, coeff)

coeff = rep(0,3)
coeff[1] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(Factor1), data=comb))$coefficients[5,1]
coeff[2] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(Factor2), data=comb))$coefficients[5,1]
coeff[3] = summary(lm(alpha_n ~ scale(iq) + scale(age) + gender + scale(Factor3), data=comb))$coefficients[5,1]
coeff = coeff*100
fact_n = data.frame(factors, coeff)


## Analysis - [alpha_p - alpha_n] / [alpha_p + alpha_n]

summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(sds_total), data=comb))
summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(oci_total), data=comb))
summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(stai_total), data=comb))
summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(lsas_total), data=comb))
summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(bis_total), data=comb))
summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(scz_total), data=comb))
summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(aes_total), data=comb))
summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(eat_total), data=comb))
summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(audit_total), data=comb))

summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(Factor1), data=comb))
summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(Factor2), data=comb))
summary(lm(((alpha_p-alpha_n)/(alpha_p+alpha_n)) ~ scale(iq) + scale(age) + gender + scale(Factor3), data=comb))


### Plotting ################

ggplot(data=symp_p, aes(x=symptoms, y=coeff, fill=symptoms)) +
  geom_bar(colour="black", stat="identity") +
  xlab("Clinical Scales") + ylab("Percent Change in Positive Learning") +
  ggtitle("Percentage Change in Positive Learning By Clinical Scale") +
  guides(fill=FALSE)

ggplot(data=fact_p, aes(x=factors, y=coeff, fill=factors)) +
  geom_bar(colour="black", stat="identity") +
  xlab("Clinical Scales") + ylab("Percent Change in Positive Learning") +
  ggtitle("Percentage Change in Positive Learning By Factor") +
  guides(fill=FALSE) 
  
ggplot(data=symp_n, aes(x=symptoms, y=coeff, fill=symptoms)) +
  geom_bar(colour="black", stat="identity") +
  xlab("Clinical Scales") + ylab("Percent Change in Negative Learning") +
  ggtitle("Percentage Change in Negative Learning By Clinical Scale") +
  guides(fill=FALSE)

ggplot(data=fact_n, aes(x=factors, y=coeff, fill=factors)) +
  geom_bar(colour="black", stat="identity") +
  xlab("Clinical Scales") + ylab("Percent Change in Negative Learning") +
  ggtitle("Percentage Change in Negative Learning By Factor") +
  guides(fill=FALSE) 
