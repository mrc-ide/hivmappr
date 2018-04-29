#' This script generates estiamtes for district-level prevalence and ART coverage from Malawi DHS and quarterly reporting data
#' 
devtools::load_all("~/Documents/Research/spatial-estimates/ART-coverage/manuscript/prevartcov/")
library(rstan)

library(maptools)

#' Load Malawi data
data(mwsh)

#' Prepare data inputs for model
standat <- with(mwsh@data, list(N_reg = length(idx),
                                district = district,
                                prev_est = prev_survey,
                                prev_se  = prev_survey_se,
                                anc1_obs = cbind(ancrt_n - ancrt_pos,
                                                 ancrt_pos - ancrt_art,
                                                 ancrt_art),
                                pop15pl_i = pop15pl,
                                pop15to49_i = pop15to49,
                                art15pl_i = adultart,
                                prev_ratio = 1.06,
                                loo_i = 0))

#' Fit model
mod <- sampling(prevartcov:::stanmodels$model, data=c(standat, model_id=5), iter=8000)

#' Extract estimates for prevalence and ART coverage

mwshest <- mwsh
mwshest$prev_mod5 <- colMeans(extract(mod, "prev_i")$prev_i)
mwshest$artcov_mod5 <- colMeans(extract(mod, "artcov_i")$artcov_i)

#' Save data in package
devtools::use_data(mwshest, overwrite=TRUE)
