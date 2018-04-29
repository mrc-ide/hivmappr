data {
  int N_reg;

  // survey
  vector[N_reg] prev_est;
  vector[N_reg] prev_se;

  int anc1_obs[N_reg, 3];

  int P_i[N_reg];
  int R_i[N_reg];
  vector[N_reg] pop15pl_i;
  vector[N_reg] pop15to49_i;
  vector[N_reg] art15pl_i;
  real prev_ratio;  // ratio of age 15+ prevalence to age 15-49 prevalence
  //
  real omega;

  // national incidence
  real log_lambda_nat_mean;
  real<lower=0> log_lambda_nat_sd;

  // prior arguments
  real<lower=0> sigma_u_sd;
}
transformed data {
  vector[N_reg] l_prev_est;
  vector[N_reg] l_prev_se;
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
}
transformed parameters{

  vector[N_reg] rho_i;
  vector[N_reg] l_rho_i;
  vector[N_reg] rho_anc_i;
  vector<lower=0, upper=1>[N_reg] alpha_i;
  vector[N_reg] l_alpha_i;
  vector<lower=0, upper=1>[N_reg] alpha_anc_i;

  rho_i = prop_art_i + (1 - prop_art_i) .* inv_logit(l_rho_i_raw); // <lower=prop_art_i, upper=1>
  l_rho_i = logit(rho_i);
  rho_anc_i = inv_logit(l_rho_i + l_rho_ancbias_i);
  
  alpha_i = prop_art_i ./ rho_i;
  l_alpha_i = logit(alpha_i);
  alpha_anc_i = inv_logit(l_alpha_i + l_alpha_ancbias_i);
}
model {

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
  l_prev_est ~ normal(l_rho_i, l_prev_se);

  // ART data likelihood
  target += log1p_exp(-l_rho_i) - log(alpha_i) - log1m(alpha_i); // propart -> l_alpha
  l_alpha_i  ~ normal(l_alpha0, sigma_l_alpha);

  // ANC data likelihood
  for(i in 1:N_reg){
    vector[3] anc1dist;
    
    anc1dist[1] = 1.0 - rho_anc_i[i];              // HIV negative
    anc1dist[2] = rho_anc_i[i] * (1 - alpha_anc_i[i]);  // HIV positive, not on ART
    anc1dist[3] = rho_anc_i[i] * alpha_anc_i[i];      // HIV positive, already on ART
    
    anc1_obs[i] ~ multinomial(anc1dist);
  }
}
generated quantities {

  // incidence
  real<lower=0> sigma_u;
  real log_kappa0;
  vector[N_reg] u_i;
  vector[N_reg] lambda_i;
  vector[N_reg]  infections_i;
  real infections;
  real rho;
  real alpha; 
  real lambda;

  sigma_u = sigma_u_sd * inv_Phi(uniform_rng(0.5, 1));
  for(i in 1:N_reg)
    u_i[i] = normal_rng(0, sigma_u);

  {
    vector[N_reg] w_i;
    vector[N_reg] S_i;
    w_i = rho_i .* (1 - omega*alpha_i) .* exp(u_i);
    S_i = (1 - rho_i) .* pop15to49_i;
    log_kappa0 = normal_rng(log_lambda_nat_mean - log(dot_product(w_i, S_i)) + log(sum(S_i)), log_lambda_nat_sd);

    lambda_i = exp(log_kappa0) * w_i;
  }

  infections_i = lambda_i .* (1 - rho_i) .* pop15to49_i;
  infections = sum(infections_i);
  rho = dot_product(rho_i, pop15to49_i) / sum(pop15to49_i);
  alpha = dot_product(alpha_i, rho_i .* pop15to49_i) / dot_product(rho_i, pop15to49_i);
  lambda = infections / dot_product(1 - rho_i, pop15to49_i);
}
