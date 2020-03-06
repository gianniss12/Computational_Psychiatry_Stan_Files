library(rstan)
library(parallel)
library(matlab)
library(rstanmulticore)

mc.cores = parallel::detectCores()

#install.packages('StanHeaders')

setwd('C:/Users/jgior/Desktop/Internship/Modeling-eLife')

data = read.csv('data/modelData_s1_forJonny.dat',sep=" ")

##
data = subset(data, trial_num<=400)

subs = unique(data$subj)

NS = length(subs)
MT=200
NT = matrix(0,NS)
c1 = array(0,dim=c(NS,MT))
c2 = array(0,dim=c(NS,MT))
st = array(0,dim=c(NS,MT))
r = array(0,dim=c(NS,MT))

    
for (i in 1:NS) {
    
     NT[i] = nrow(subset(data,subj==subs[i]));
     c1[i,1:NT[i]] = subset(data,subj==subs[i])$step1resp - 1;
     c2[i,1:NT[i]] = subset(data,subj==subs[i])$step2resp - 1;
     st[i,1:NT[i]] = subset(data,subj==subs[i])$step2state - 1;
     r[i,1:NT[i]] = subset(data,subj==subs[i])$reward;
    
}

standata = list(NS= NS, MT=MT, NT= as.vector(NT), c1=c1, c2=c2, st=st, r=r)

seed <- runif(1,1,1e6); # Stan wants a random integer from 1 to max supportable

fit <- stan(file = 'stan_fitting_scripts/JG_stan_sepAlphas_sepStages_noLambda_noDecay_noRescaling.stan', data = standata, seed=seed, iter = 1, chains = 1)

# sflist = mclapply(1:4, mc.cores = 4, function(i) stan(fit=fit, seed = sample.int(.Machine$integer.max, 1), data = standata, iter=1000, chains = 1, chain_id = i)) #, refresh = -1))

# Mac or Linux
#sflist = mclapply(1:4, mc.cores = 1, function(i) stan(fit=fit, seed = sample.int(.Machine$integer.max, 1), data = standata, iter=4000, chains = 1, chain_id = i , refresh = 1))
#fit <- sflist2stanfit(sflist)


# Windows
fit = pstan(fit=fit, seed = sample.int(.Machine$integer.max, 1), data = standata, iter=4000, chains = 4, refresh = 1)

save.image(file = "stan_fits/JG_s1_Sharp_noDecay_March4_2020.RData", version = NULL, ascii = FALSE)


# run this part first and then save how you wish later
if(FALSE){
m = extract(fit, permuted=TRUE) #, permuted = FALSE, inc_warmup = FALSE) # return an array

#
print(c('MB', median(m$b1mm), quantile(m$b1mm, c(0.025, 0.975) )))

bySubj = data.frame(subj=unique(data$subj), beta1m=rep(0,(length(unique(data$subj)))), betac=rep(0,(length(unique(data$subj)))), lambda=rep(0,(length(unique(data$subj)))), alpha1=rep(0,(length(unique(data$subj)))), rHat=rep(0,(length(unique(data$subj)))),w=rep(0,(length(unique(data$subj)))), beta1t=rep(0,(length(unique(data$subj)))))

for( i in 1:90){
    bySubj[i,2:7] = c(mean(m$beta1m[,i]),mean(m$betac[,i]),mean(m$lambda[,i]),mean(m$alpha1[,i]), mean(m$rhat[,i]), mean(m$w[,i]), mean(m$beta1t[,i]))
    
}


# m1 is a model with: dev_diff_all_z + decays Q-values for non-chosen + stickiness is common across both states

m1 <- fit
merge_m1 = merge(merged, bySubj, by = "subj", all=TRUE)
}


