###############################################################
###############################################################
#### Running Stan Model for Two Step Computational Model ######
###############################################################
###############################################################
# Seow 2017
# Runs Stan Model for Subject level betas
# Extract the betas (both separate MF/MB or w)
# Compares betas with regression coeffients


rm(list=ls())

library(rstan)
library(parallel)
library(matlab)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

dataTemp = read.csv('TS_CompModelPrep2.csv',sep=",",header=T)

colnames(dataTemp) <- c("subject","trial_num","step2resp","step1resp","step2state","reward")
#dataTemp$rewardCh <- ifelse(dataTemp$reward==1, 1,-1 )
dataTemp <- dataTemp[order(dataTemp$subj,  dataTemp$trial_num), ]
dataTemp <- droplevels(dataTemp)

#remove na trials
data =dataTemp[complete.cases(dataTemp), ]

#remove here #"EEGTDPRJ0104F391" "EEGTDPRJ31F6E6A6" "EEGTDPRJ7AA3D990" "EEGTDPRJDE9CB11A" to match eeg data
data<-data[!(data$subj=="EEGTDPRJ0104F391"| data$subj== "EEGTDPRJ31F6E6A6"| data$subj== "EEGTDPRJ7AA3D990"| data$subj== "EEGTDPRJDE9CB11A"),]
data$subj<-droplevels(data$subj)

set_cppo("fast")
data<- data[order(data[,1], data[, 2]),] #make sure trial num is ordered, just in case

subs = unique(data$subject)

NS = length(subs)
MT = 150
NT = matrix(0,NS)
c1 = array(0,dim=c(NS,MT))
c2 = array(0,dim=c(NS,MT))
st = array(0,dim=c(NS,MT))
r = array(0,dim=c(NS,MT))


for (i in 1:NS) {
  
  NT[i] = nrow(subset(data,subject==subs[i]));
  c1[i,1:NT[i]] = subset(data,subject==subs[i])$step1resp - 1; # 1 or 0
  c2[i,1:NT[i]] = subset(data,subject==subs[i])$step2resp - 1; # 1 or 0
  st[i,1:NT[i]] = subset(data,subject==subs[i])$step2state - 1; # 1 or 2
  r[i,1:NT[i]] = subset(data,subject==subs[i])$reward; # 1 or 0
  
}


standata = list(NS= NS, MT=MT, NT= as.vector(NT), c1=c1, c2=c2, st=st, r=r)

seed <- runif(1,1,1e6); # Stan wants a random integer from 1 to max supportable

#fit <- stan(file = 'TS_stan_withdecay_sepBetas.stan', data = standata, seed=seed, iter = 1, chains = 1, control = list(adapt_delta = 0.99))
#fit = stan(file='TS_stan_withdecay_sepBetas.stan', data=standata,verbose=FALSE,save_warmup=FALSE, seed = sample.int(.Machine$integer.max, 1),iter=4000,control=list(adapt_delta=0.99,stepsize=.01))
#fit = stan(file='TS_stan_withdecay_sepBetasWLambda3.stan', data=standata,verbose=FALSE,save_warmup=FALSE, seed = sample.int(.Machine$integer.max, 1),iter=4000,control=list(adapt_delta=0.99,stepsize=.01))

#fit = stan(file='TS_stan_withdecay_sepBetasWLambda3.stan', data=standata,verbose=FALSE,save_warmup=FALSE, seed = sample.int(.Machine$integer.max, 1),iter=4000,chains = 4, cores = 4)

# Mac OS update to Catalina doesn't let me run multiple cores anymore..
fit <- stan(file = 'TS_stan_withdecay_sepBetasWLambda2-2.stan', data = standata, seed=seed, iter = 1, chains = 1)
sflist = mclapply(1:4, mc.cores = 4, function(i) stan(fit=fit, seed = sample.int(.Machine$integer.max, 1), data = standata, iter=4000, chains = 1, chain_id = i)) #, refresh = -1))

#iters should be 4000
fit <- sflist2stanfit(sflist)

save.image(file = "TS_stan_withdecay_sepBetasWLambda2-2.RData", version = NULL, ascii = FALSE)


######################################################
#### Extracting Betas #################################
# #####################################################
library("bayesplot")
library("ggplot2")
library("rstan")
library(parallel)
library(matlab)
library(shinystan)


# Extract fit object
m <- extract(fit, permuted=TRUE)
str(m)
code <- get_stancode(fit) # the stan code
cat(code)

# Check beta effect samples and rhat (summary)
beta1MF_summary=summary(fit,pars=c('beta1t'),probs=c(0.5))$summary
beta1MB_summary=summary(fit,pars=c('beta1m'),probs=c(0.5))$summary
pers_summary=summary(fit,pars=c('betac'),probs=c(0.5))$summary


## Some MCMC Diagnositics
lp_cp <- log_posterior(fit)
head(lp_cp)
np_cp <- nuts_params(fit)
head(np_cp)

# Plot rHats
rhats <- rhat(fit, pars = "beta1m")
print(rhats)
color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(rhats)  #+ yaxis_text(hjust = 1)


# Plot effSam
ratios_cp <- neff_ratio(fit, pars = "beta1m")
print(ratios_cp)
mcmc_neff(ratios_cp, size = 2)


# Plot invividual trace plots to check convergence
#traceplot(fit,pars ='beta1m[1]')

# creating empty array for betas (if all is great)
bySubj = data.frame(subj=unique(data$subj), betac=rep(0,,(length(unique(data$subj)))), alpha1=rep(0,,(length(unique(data$subj)))), beta1m=rep(0,,(length(unique(data$subj)))), beta1t=rep(0,,(length(unique(data$subj)))), beta2=rep(0,,(length(unique(data$subj)))) , lambda=rep(0,,(length(unique(data$subj))))   )

 for( i in 1:length(unique(data$subj))){
   bySubj[i,2:7] = c(mean(m$betac[,i]),mean(m$alpha1[,i]), mean(m$beta1m[,i]),mean(m$beta1t[,i]), mean(m$beta2[,i]), mean(m$lambda[,i]))
 }

write.csv(bySubj, file = "TS_stan_withdecay_sepBetasWLambdaTestR.csv",row.names=FALSE)
