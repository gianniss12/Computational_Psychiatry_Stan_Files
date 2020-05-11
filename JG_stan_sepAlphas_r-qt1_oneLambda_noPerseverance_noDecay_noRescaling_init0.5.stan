// Title: JG_stan_sepAlphas_r-qt1_oneLambda_noPerseverance_noDecay_noRescaling
// By: Jonny Giordano
// Date: April 27th, 2020

// Purpose: This script fits multiple learning rates for positive, negative valence
// and is based on the M3 version by Brown, but does not include the stickiness parameter

// Notes:
// This is a simplified version of the Daw model. No seperate parameters for each stage
// 1) assume lambda = 1
// 2) assume multiple learning rates
// 3) code separate model based beta and model free betaa
// 4) DO NOT rescale the learning rule updates (by alpha) to get more consistent scaling on betas
// 5) DO NOT decay unchosen options
// 6) DO NOT include Stickiness (aka Perseverance)


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
  real<lower=0,upper=1> apm;
  real<lower=0> apss;
  
  real<lower=0,upper=1> anm;
  real<lower=0> anss;

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
  # real bcm;
  # real<lower=0> bcs;

  // per subject learning rates
  //real<lower=0,upper=1> alpha1[NS];
  real<lower=0,upper=1> alpha_p[NS];
  real<lower=0,upper=1> alpha_n[NS];

  // per subject betas
  real beta1m[NS];
  real beta1t[NS];
  real beta2[NS];
  # real betac[NS];    

 
}
transformed parameters {
	real apa;
	real apb;
	
	real ana;
	real anb;
	
    
  //for learning rate
	apa = apm * pow(apss,-2);
	apb = pow(apss,-2) - apa;
	
	ana = anm * pow(anss,-2);
	anb = pow(anss,-2) - ana;
	
   
}
model {

  // priors
  
  #bcm ~ normal(0,100);
  #bcs ~ cauchy(0,2.5); //cauchy(0,2.5);cauchy(0,2.5);// above this was .1 in w version...? not sure why.
  
  b1mm ~ normal(0,100);
  b1ms ~ cauchy(0,2.5);
  b1tm ~ normal(0,100);
  b1ts ~ cauchy(0,2.5);
  
  b2m ~ normal(0,100);
  b2s ~ cauchy(0,2.5);

  // group level model (parameters per subject)

	for (s in 1:NS) {
        #int pc;
        int tcounts[2,2];
        real qm[2];
        real qt1[2];
        real qt2[2,2];
        real delta_1;
        real delta_2;
        real delta_diff;
  
      alpha_p[s] ~ beta(apa,apb);
      alpha_n[s] ~ beta(ana, anb);
      beta1m[s] ~ normal(b1mm,b1ms);
      beta1t[s] ~ normal(b1tm,b1ts);
      beta2[s] ~ normal(b2m,b2s);
      # betac[s] ~ normal(bcm,bcs); Remove "stickiness"" beta


      tcounts = rep_array(0, 2, 2); //// square matrix: choice x 2nd Stage State, Reset to all zeros
      qm = rep_array(0.5, 2); //// Qm = Model Based Value, Rest to zeros
      qt1 = rep_array(0.5, 2); //// Qt Model Free Value, First Stage , Reset to zeros
      qt2 = rep_array(0.5, 2, 2); //// Model Free Values, Square Matrix:  2nd Stage State x Final Choice, Reset to zeros
      #pc = 0;
      

      for (t in 1:NT[s]) {
        int nc1;
        int nc2;
        int nst;

        qm[1] = (int_step(tcounts[1,1]+tcounts[2,2]-tcounts[1,2]-tcounts[2,1]) ? fmax(qt2[1,1],qt2[1,2]) : fmax(qt2[2,1],qt2[2,2]));
        
        qm[2] = (int_step(tcounts[1,1]+tcounts[2,2]-tcounts[1,2]-tcounts[2,1]) ? fmax(qt2[2,1],qt2[2,2]) : fmax(qt2[1,1],qt2[1,2]));
          
          
        c1[s,t] ~ bernoulli_logit((beta1m[s])  * (qm[2] - qm[1])
          + (beta1t[s])  * (qt1[2] - qt1[1]));
          # + (betac[s]) * pc ); Remove stickiness beta
        
        c2[s,t] ~ bernoulli_logit((beta2[s]) * (qt2[st[s,t],2] - qt2[st[s,t],1]));
        
        ## Comment out so perseverance will remain 0
        #pc = 2 * c1[s,t] - 1;
        
        tcounts[c1[s,t]+1,st[s,t]] = tcounts[c1[s,t]+1,st[s,t]] + 1;
        
        // delta_1 = qt2 - qt1
        delta_1 = qt2[st[s,t],c2[s,t]+1] - qt1[c1[s,t]+1];
        
        // delta_2 = r - qt2
        delta_2 = r[s,t] - qt2[st[s,t],c2[s,t]+1];
        
        // delta_diff = delta_1 - delta2
        delta_diff = r[s,t] - qt1[c1[s,t]+1];
        
        // Multiple Learning Rates : qt1 = qt1 + alpha_x * (r-qt1)
        qt1[c1[s,t]+1] = (delta_diff >= 0) ? qt1[c1[s,t]+1] + (alpha_p[s] * delta_diff) : qt1[c1[s,t]+1] + (alpha_n[s] * delta_diff);

        // Multiple Learning Rates : qt2 = qt2 + alpha_x * (r-qt2)
        qt2[st[s,t],c2[s,t]+1] = (delta_2 >= 0) ? qt2[st[s,t],c2[s,t]+1] + (alpha_p[s] * delta_2) : qt2[st[s,t],c2[s,t]+1] + (alpha_n[s] * delta_2);


      }
      
    }

}
