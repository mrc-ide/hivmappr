data {

  int N_reg;

  vector[N_reg] pop15pl_i;
  vector[N_reg] pop15to49_i;
  vector[N_reg] art15pl_i;
  real prev_ratio;  // ratio of age 15+ prevalence to age 15-49 prevalence

  // survey prevalence
  int Nobs_prev;               // Number of observations
  int idx_prev[Nobs_prev];     // region index
  vector[Nobs_prev] prev_est;
  vector[Nobs_prev] prev_se;

  // routine ANC testing
  int anc1_obs[N_reg, 3];

  // recent infection testing
  int Nobs_rec;
  int idx_rec[Nobs_rec];
  int P_i[Nobs_rec];
  int R_i[Nobs_rec];

  //
  real OmegaT0;
  real sigma_OmegaT;
  real betaT0;
  real sigma_betaT;
  real T;
  //
  real omega;

  // covariate for incidence
  int Nkappa;
  matrix[N_reg, Nkappa] Xkappa;

  // prior arguments
  real<lower=0> sigma_u_sd;

}
transformed data {
  vector[Nobs_prev] l_prev_est;
  vector[Nobs_prev] l_prev_se;
  vector[N_reg] prop_art_i; // proportion on ART, lower bound for prevalence based on constraint prev_i * prev_ratio > art_i / pop15pl_i;

  l_prev_est = logit(prev_est);
  l_prev_se = prev_se ./ (prev_est .* (1 - prev_est));
  prop_art_i = art15pl_i ./ (prev_ratio * pop15pl_i);
}
parameters {

  // prevalence
  real l_rho0;
  real<lower=0> sigma_l_rho;
  vector[N_reg] l_rho_i_raw;
  
  real l_rho_ancbias0;
  real<lower=0> sigma_rho_ancbias;
  vector[N_reg] l_rho_ancbias_i;

  // ART coverage
  real l_alpha0;
  real<lower=0> sigma_l_alpha;
  
  vector[N_reg] l_alpha_ancbias_i;
  real l_alpha_ancbias0;
  real<lower=0> sigma_alpha_ancbias;

  // incidence
  vector[Nkappa] kappa;
  real<lower=0> sigma_u;
  vector[N_reg] u_raw;
  real OmegaT;
  real<lower=0> betaT_raw;
}
transformed parameters{

  vector[N_reg] rho_i;
  vector[N_reg] l_rho_i;
  vector[N_reg] rho_anc_i;
  vector<lower=0, upper=1>[N_reg] alpha_i;
  vector[N_reg] l_alpha_i;
  vector<lower=0, upper=1>[N_reg] alpha_anc_i;

  real betaT;
  vector[N_reg] u_i;
  
  betaT = betaT_raw * sigma_betaT + betaT0;
  u_i = sigma_u * u_raw;

  rho_i = prop_art_i + (1 - prop_art_i) .* inv_logit(l_rho_i_raw); // <lower=prop_art_i, upper=1>
  l_rho_i = logit(rho_i);
  rho_anc_i = inv_logit(l_rho_i + l_rho_ancbias_i);
  
  alpha_i = prop_art_i ./ rho_i;
  l_alpha_i = logit(alpha_i);
  alpha_anc_i = inv_logit(l_alpha_i + l_alpha_ancbias_i);
}
model {

  vector[N_reg] kappa_i;
  vector[N_reg] pR;


  // priors

  l_rho0 ~ normal(-2, 5);
  sigma_l_rho ~ normal(0, 2.5);

  target += l_rho_i_raw - log1p_exp(l_rho_i_raw) - log(rho_i); // rho_i_raw -> l_rho_i
  l_rho_i ~ normal(l_rho0, sigma_l_rho);

  l_alpha0 ~ normal(0, 5);
  sigma_l_alpha ~ normal(0, 2.5);

  l_rho_ancbias0 ~ normal(0, 5);
  sigma_rho_ancbias ~ normal(0, 2.5);
  l_rho_ancbias_i ~ normal(l_rho_ancbias0, sigma_rho_ancbias);

  l_alpha_ancbias0 ~ normal(0, 5);
  sigma_alpha_ancbias ~ normal(0, 2.5);
  l_alpha_ancbias_i ~ normal(l_alpha_ancbias0, sigma_alpha_ancbias);
  

  // prevalence likelihood
  l_prev_est ~ normal(l_rho_i[idx_prev], l_prev_se);

  // ART data likelihood
  target += log(prop_art_i) - log1m(alpha_i); // propart -> l_alpha
  l_alpha_i  ~ normal(l_alpha0, sigma_l_alpha);

  // ANC data likelihood
  for(i in 1:N_reg){
    vector[3] anc1dist;
    
    anc1dist[1] = 1.0 - rho_anc_i[i];              // HIV negative
    anc1dist[2] = rho_anc_i[i] * (1 - alpha_anc_i[i]);  // HIV positive, not on ART
    anc1dist[3] = rho_anc_i[i] * alpha_anc_i[i];      // HIV positive, already on ART
    
    anc1_obs[i] ~ multinomial(anc1dist);
  }

  // recent infection
  
  kappa ~ normal(0, 1);
  
  u_raw ~ normal(0, 1); // u_i ~ normal(0, sigma_u)
  sigma_u ~ normal(0, sigma_u_sd); // half-normal(0, sigma_u_sd) prior: constrain
  OmegaT ~ normal(OmegaT0, sigma_OmegaT);
  betaT_raw  ~ normal(0, 1); // betaT ~ normal(betaT0, sigma_betaT)

  kappa_i = exp(Xkappa * kappa + log1m(omega*alpha_i) + u_i);
  pR =  kappa_i .* (1-rho_i) * (OmegaT - betaT*T) + betaT;
  R_i ~ binomial(P_i, pR[idx_rec]);
}
generated quantities {
  vector[N_reg] lambda_i;
  vector[N_reg]  infections_i;
  vector[N_reg] excess_i;
  real infections;
  real rho;
  real alpha; 
  real lambda;
  //
  excess_i = Xkappa * kappa + u_i - mean(Xkappa * kappa + u_i);
  {
    vector[N_reg] kappa_i;
    kappa_i = exp(Xkappa * kappa + log1m(omega*alpha_i) + u_i);
    lambda_i = kappa_i .* rho_i;
  }
  infections_i = lambda_i .* (1 - rho_i) .* pop15to49_i;
  infections = sum(infections_i);
  rho = dot_product(rho_i, pop15to49_i) / sum(pop15to49_i);
  alpha = dot_product(alpha_i, rho_i .* pop15to49_i) / dot_product(rho_i, pop15to49_i);
  lambda = infections / dot_product(1 - rho_i, pop15to49_i);
}
