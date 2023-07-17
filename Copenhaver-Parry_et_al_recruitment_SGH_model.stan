////Negative Binomial Regression Model
////Copenhaver-Parry, P.E., Byerly, S.N., Greenwood, J.W., Retter, C.J., Beedlow, P.A., & Lee, E.H.
////Simultaneous changes in climate and competition limit forest recruitmetn with increasing climatic stress

data{
	int<lower=1> s;					//number of sites
	int<lower=1> t;                 		//number of years
	int<lower=0> N;					//number of seedling count observations
	int<lower=1, upper=s> site[N];  		//vector of site numbers (one value for each observation)
	int<lower=1, upper=t> year[N];			//vector of year numbers (one value for each observation)
	int<lower=0> y[N];				//response variable (seedling counts)
	int<lower=1> K;                 		//number of covariates
	matrix[N,K] X;                  		//covariate matrix
}

parameters{
  vector[s] alpha_s;              			//random site effect
  vector[t] alpha_t;              			//random year effect
  vector[K] beta[s];               			//matrix of site-varying regression coefficients; row for each site, column for each coefficient
  real mu_s;                      			//mean of random site effects
  real<lower=0> tau_s;                     		//sd of random site effects
  real mu_t;                      			//mean of random year effects
  real<lower=0> tau_t;                     		//sd of random year effects
  vector[K] beta_mu;                   			//vector of means of site-varying covariate effects
  vector<lower=0>[K] tau;                   		//vector of sds of site-varying covariate effects
  real<lower=0> reciprocal_r;
}

transformed parameters{
  real log_p[N];
  real r;
  for (n in 1:N){
	  log_p[n]=alpha_s[site[n]] + alpha_t[year[n]] + X[n]*beta[site[n]];
	}
 r = 1. / reciprocal_r;
}

model{
  vector[N]mu;                 
  //priors
  alpha_s~normal(mu_s,tau_s);
  alpha_t~normal(mu_t,tau_t);
  beta_mu~normal(0,2);
  for(S in 1:s){
    	beta[S]~normal(beta_mu,tau);
  }
  tau~gamma(2,0.1);
  mu_s~normal(0,3);
  mu_t~normal(0,3);
  tau_s~gamma(2,0.1);
  tau_t~gamma(2,0.1);
  reciprocal_r~cauchy(0.,5);
  //likelihood
  for(n in 1:N){
	target +=neg_binomial_2_log_lpmf(y[n]|log_p[n], r);
  }
}

generated quantities {
	vector[N] y_pred;				//predicted counts
	vector[N] log_lik;				//log-likelihood for model evaluation
	 for(n in 1:N){								
    	y_pred[n]=neg_binomial_2_log_rng(log_p[n],r);			        
    	log_lik[n]=neg_binomial_2_log_lpmf(y[n]|log_p[n],r);	
  }
}
