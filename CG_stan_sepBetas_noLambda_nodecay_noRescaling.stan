# this is a stan version of the M3 version used by Vanessa Brown
# 1) assume lambda = 1
# 2) assume a single learning rate (alpha)
# 3) code separate model based beta and model free betaa
# 4) DO NOT rescale the learning rule updates (by alpha) to get more consistent scaling on betas
# 5) DO NOT decay unchosen options


data {
	int NS;  // number of subjects
	int MT;  // maximum trials per subject
    
  int r[NS,MT];   // 1 / 0 (changed my stan.r file to make this 1,-1 as required)
  int NT[NS];  // actual numbers of trials per subject

	int c1[NS,MT];  // 0 / 1
	int c2[NS,MT];  // 0 / 1
	int st[NS,MT];  // 1 / 2

}

# this lays out all the model parameters -- both the group level and individual subject level

parameters {
 
  # parameterize hyperpriors on this distribution for LR not in terms of a, b directly, but instead 
  # mean (a/(a+b)) and a function of sample size 1/sqrt(a+b), which is roughly related to its standard deviation
  # the (implicit, improper) infinite uniform prior on these parameters is then appropriate
  # see Gelman et al p. 128 & exercise 5.7

  # group level learning rate (no constraint here, that is done later)
  #real a1m; // no constraint
  real<lower=0,upper=1> a1m;
  real<lower=0> a1ss;

  # group level model-based beta (mean and SD)
  real b1mm;
  real<lower=0> b1ms;

  # group level TD-1 beta 
  real b1tm;
  real<lower=0> b1ts;

  # group level beta for second stage choices
  real b2m;
  real<lower=0> b2s;

  # group level perseveration beta
  real bcm;
  real<lower=0> bcs;

  # per subject learning rates
  real<lower=0,upper=1> alpha1[NS];

  # per subject betas
  real beta1m[NS];
  real beta1t[NS];
  real beta2[NS];
  real betac[NS];    

 
}
transformed parameters {
	real a1a;
	real a1b;	
    
  #for learning rate
	a1a <- a1m * pow(a1ss,-2);
	a1b <- pow(a1ss,-2) - a1a;
	
   
}
model {

  # priors
  
  bcm ~ normal(0,100);
  bcs ~ cauchy(0,2.5); #cauchy(0,2.5);cauchy(0,2.5);# above this was .1 in w version...? not sure why.
  
  b1mm ~ normal(0,100);
  b1ms ~ cauchy(0,2.5);
  b1tm ~ normal(0,100);
  b1ts ~ cauchy(0,2.5);
  
  b2m ~ normal(0,100);
  b2s ~ cauchy(0,2.5);

  # group level model (parameters per subject)

	for (s in 1:NS) {
        int pc;
        int tcounts[2,2];
        real qm[2];
        real qt1[2];
        real qt2[2,2];
        real delta_1;
        real delta_2;
  
      alpha1[s] ~ beta(a1a,a1b);
      beta1m[s] ~ normal(b1mm,b1ms);
      beta1t[s] ~ normal(b1tm,b1ts);
      beta2[s] ~ normal(b2m,b2s);
      betac[s] ~ normal(bcm,bcs);

      for (i in 1:2) for (j in 1:2) tcounts[i,j] <- 0;
      for (i in 1:2) {qm[i] <- 0; qt1[i] <- 0;}
      for (i in 1:2) for (j in 1:2) qt2[i,j] <- 0;
      pc <- 0;

      for (t in 1:NT[s]) {
        int nc1;
        int nc2;
        int nst;

        qm[1] <- if_else(int_step(tcounts[1,1]+tcounts[2,2]-tcounts[1,2]-tcounts[2,1]), 
          fmax(qt2[1,1],qt2[1,2]), fmax(qt2[2,1],qt2[2,2]));

        qm[2] <- if_else(int_step(tcounts[1,1]+tcounts[2,2]-tcounts[1,2]-tcounts[2,1]), 
          fmax(qt2[2,1],qt2[2,2]), fmax(qt2[1,1],qt2[1,2]));
          #previously betac[s]
          
        c1[s,t] ~ bernoulli_logit((beta1m[s])  * (qm[2] - qm[1])
          + (beta1t[s])  * (qt1[2] - qt1[1])
          + (betac[s]) * pc );
        
        c2[s,t] ~ bernoulli_logit((beta2[s]) * (qt2[st[s,t],2] - qt2[st[s,t],1]));

        pc <- 2 * c1[s,t] - 1;
        tcounts[c1[s,t]+1,st[s,t]] <- tcounts[c1[s,t]+1,st[s,t]] + 1;
        
        delta_1 = qt2[st[s,t],c2[s,t]+1]-qt1[c1[s,t]+1]; 
        delta_2 = r[s,t] - qt2[st[s,t],c2[s,t]+1];
        
        ## first step updating TD value: (old qt1 value + (value of second stage state - value of choice) * alpha)
        ## update first state: (old qt1 value + (reward - value of second state)*alpha)
        qt1[c1[s,t]+1] = qt1[c1[s,t]+1] + alpha1[s]*(delta_1+delta_2);
        
        ## second step updating TD value: (old qt2 value + (reward - value of second state)* alpha)
        qt2[st[s,t],c2[s,t]+1] <- qt2[st[s,t],c2[s,t]+1] + (r[s,t] - qt2[st[s,t],c2[s,t]+1] ) * alpha1[s];
        
        
      }
      
    }

}
