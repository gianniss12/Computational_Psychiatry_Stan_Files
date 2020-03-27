// Title: stan_sepAlphas_sepStages_noLambda_noDecay_noRescaling
// By: Jonny Giordano
// Data: February 26th, 2020

// Purpose: This script fits multiple learning rates for positive, negative valence at each stage
// and is based on the M3 version by Brown

// Notes:
// This is a simplified version of the Daw model. No seperate parameters for each stage
// 1) assume lambda = 1
// 2) assume multiple learning rates
// 3) code separate model based beta and model free betaa
// 4) DO NOT rescale the learning rule updates (by alpha) to get more consistent scaling on betas
// 5) DO NOT decay unchosen options


data {
	int NS;  //// number of subjects
	int MT;  //// maximum trials per subject
    
  int r[NS,MT];   //// 1 // 0 (changed my stan.r file to make this 1,-1 as required)
  int NT[NS];  //// actual numbers of trials per subject

	int c1[NS,MT];  //// 0 // 1
	int c2[NS,MT];  //// 0 // 1
	int st[NS,MT];  //// 1 // 2

}

// this lays out all the model parameters -- both the group level and individual subject level

parameters {
 
  // parameterize hyperpriors on this distribution for LR not in terms of a, b directly, but instead 
  // mean (a//(a+b)) and a function of sample size 1//sqrt(a+b), which is roughly related to its standard deviation
  // the (implicit, improper) infinite uniform prior on these parameters is then appropriate
  // see Gelman et al p. 128 & exercise 5.7

  // group level learning rate (no constraint here, that is done later)
  //real a1m; //// no constraint\
  // Positive, Negative Alphas
  real<lower=0,upper=1> a1pm; // Positive Alpha, Stage 1
  real<lower=0> a1pss;
  
  real<lower=0,upper=1> a2pm; // Positive Alpha, Stage 2
  real<lower=0> a2pss;
  
  real<lower=0,upper=1> a1nm; // Negative Alpha, Stage 1
  real<lower=0> a1nss;
  
  real<lower=0,upper=1> a2nm; // Negative Alpha, Stage 2
  real<lower=0> a2nss;

  // group level model-based beta (mean and SD)
  real b1mm;
  real<lower=0> b1ms;

  // group level TD-1 beta 
  real b1tm;
  real<lower=0> b1ts;

  // group level beta for second stage choices
  real b2m;
  real<lower=0> b2s;

  // group level perseveration beta
  real bcm;
  real<lower=0> bcs;

  // per subject learning rates
  //real<lower=0,upper=1> alpha1[NS];
  real<lower=0,upper=1> alpha1_p[NS];
  real<lower=0,upper=1> alpha2_p[NS];
  real<lower=0,upper=1> alpha1_n[NS];
  real<lower=0,upper=1> alpha2_n[NS];

  // per subject betas
  real beta1m[NS];
  real beta1t[NS];
  real beta2[NS];
  real betac[NS];    

}

transformed parameters {
	real a1pa;
	real a1pb;
	
	real a2pa;
	real a2pb;
	
	real a1na;
	real a1nb;

	real a2na;
	real a2nb;
    
  //for learning rate
	a1pa = a1pm * pow(a1pss,-2);
	a1pb = pow(a1pss,-2) - a1pa;
	
	a2pa = a2pm * pow(a2pss,-2);
	a2pb = pow(a2pss,-2) - a2pa;
	
	a1na = a1nm * pow(a1nss,-2);
	a1nb = pow(a1nss,-2) - a1na;
	
	a2na = a2nm * pow(a2nss,-2);
	a2nb = pow(a2nss,-2) - a2na;
   
}
model {

  // priors
  
  bcm ~ normal(0,100);
  bcs ~ cauchy(0,2.5); //cauchy(0,2.5);cauchy(0,2.5);// above this was .1 in w version...? not sure why.
  
  b1mm ~ normal(0,100);
  b1ms ~ cauchy(0,2.5);
  b1tm ~ normal(0,100);
  b1ts ~ cauchy(0,2.5);
  
  b2m ~ normal(0,100);
  b2s ~ cauchy(0,2.5);

  // group level model (parameters per subject)

	for (s in 1:NS) {
        int pc;
        int tcounts[2,2];
        real qm[2];
        real qt1[2];
        real qt2[2,2];
        real delta_1;
        real delta_2;
        real delta_diff;
  
      alpha1_p[s] ~ beta(a1pa,a1pb);
      alpha2_p[s] ~ beta(a2pa,a2pb);
      alpha1_n[s] ~ beta(a1na,a1nb);
      alpha2_n[s] ~ beta(a2na,a2nb);
      beta1m[s] ~ normal(b1mm,b1ms);
      beta1t[s] ~ normal(b1tm,b1ts);
      beta2[s] ~ normal(b2m,b2s);
      betac[s] ~ normal(bcm,bcs);

      tcounts = rep_array(0, 2, 2); //// square matrix: choice x 2nd Stage State, Reset to all zeros
      qm = rep_array(0.0, 2); //// Qm = Model Based Value, Rest to zeros
      qt1 = rep_array(0.0, 2); //// Qt Model Free Value, First Stage , Reset to zeros
      qt2 = rep_array(0.0, 2, 2); //// Model Free Values, Square Matrix:  2nd Stage State x Final Choice, Reset to zeros
      pc = 0;
      

      for (t in 1:NT[s]) {
        int nc1;
        int nc2;
        int nst;

        //// Optimization
        qm[1] = (int_step(tcounts[1,1]+tcounts[2,2]-tcounts[1,2]-tcounts[2,1]) ? fmax(qt2[1,1],qt2[1,2]) : fmax(qt2[2,1],qt2[2,2]));
        
        qm[2] = (int_step(tcounts[1,1]+tcounts[2,2]-tcounts[1,2]-tcounts[2,1]) ? fmax(qt2[2,1],qt2[2,2]) : fmax(qt2[1,1],qt2[1,2]));
          
          
        c1[s,t] ~ bernoulli_logit((beta1m[s])  * (qm[2] - qm[1])
          + (beta1t[s])  * (qt1[2] - qt1[1])
          + (betac[s]) * pc );
        
        c2[s,t] ~ bernoulli_logit((beta2[s]) * (qt2[st[s,t],2] - qt2[st[s,t],1]));

        pc = 2 * c1[s,t] - 1;
        tcounts[c1[s,t]+1,st[s,t]] = tcounts[c1[s,t]+1,st[s,t]] + 1;
        
        // delta_1 = qt2 - qt1
        delta_1 = qt2[st[s,t],c2[s,t]+1] - qt1[c1[s,t]+1];
        
        // delta_2 = r - qt2
        delta_2 = r[s,t] - qt2[st[s,t],c2[s,t]+1];
        
        // delta_diff = delta_1 - delta2
        delta_diff = r[s,t] - qt1[c1[s,t]+1];
        
        // Multiple Learning Rates : qt1 = qt1 + alpha_x * (delta_1)
        qt1[c1[s,t]+1] = (delta_1 >= 0) ? qt1[c1[s,t]+1] + (alpha1_p[s] * delta_1) : qt1[c1[s,t]+1] + (alpha1_n[s] * delta_1);

        // Multiple Learning Rates : qt2 = qt2 + alpha_x * (delta_2)
        qt2[st[s,t],c2[s,t]+1] = (delta_2 >=0) ? qt2[st[s,t],c2[s,t]+1] + (alpha2_p[s] * delta_2) : qt2[st[s,t],c2[s,t]+1] + (alpha2_n[s] * delta_2);
        
        

      }
      
    }

}
