setwd("~/Documents/Research/spatial-estimates/hivincid/analysis/")
devtools::load_all()
library(rstan)
library(data.table)

rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())


## mw@data <- data.table(mw@data)

## mw@data[district %in% c("Mchinji", "Phalombe")]
           
## mwsim <- mw
## mwsim <- merge(mwsim, data.frame(district = mwsim$district,
##                                  rho_sim = summary(fit, "rho_i")$summary[,1],
##                                  alpha_sim = summary(fit, "alpha_i")$summary[,1],
##                                  lambda_sim = summary(fit, "lambda_i")$summary[,1]))

## mwsim$u_i <- ifelse(mwsim$district %in% c("Mchinji", "Phalombe"), 0.7, -0.054)
## mwsim$lambda_i <- with(mwsim@data, exp(-2.8) * rho_sim * (1-0.7*alpha_sim) * exp(u_i))
## mwsim$p_r_i <- with(mwsim@data, lambda_i * (1-rho_sim) / rho_sim * OmegaT)
## with(mwsim@data, sum(lambda_i * pop15to49 * (1-prev_mod5)) / sum(pop15to49 * (1-prev_mod5)))



## devtools::use_data(mwsim)

#' Generate simulated data

data(mwsim)
OmegaT <- 130/365

hotid <- c(11, 22)
coldid <- c(1:10, 12:21, 23:28)


data <- with(mwsim@data, list(N_reg = length(district),
                              district = district,
                              prev_est = npos / nsamp,
                              prev_se = sqrt(2 * npos * (nsamp - npos) / nsamp^3),
                              anc1_obs = cbind(ancrt_neg = ancrt_n - ancrt_pos, ancrt_noart=ancrt_pos_noart, ancrt_art=ancrt_art),
                              P_i = npos,
                              R_i = nrecent,
                              pop15pl_i = pop15pl,
                              pop15to49_i = pop15to49,
                              art15pl_i = adultart,
                              prev_ratio = 1.06,
                              OmegaT0 = 130 / 365,
                              sigma_OmegaT = ((142-118)/365)/(2*qnorm(0.975)),
                              betaT0 = 0.0,
                              sigma_betaT = 0.00001,
                              omega = 0.7,
                              T = 1.0,
                              Nkappa = 1,
                              Xkappa = matrix(10, length(district), 1)))


sim_hotspot <- function(E_recent=20){
  prev <- with(mwsim@data, sum(rho_sim * pop15to49) / sum(pop15to49))
  incid <- with(mwsim@data, sum(lambda_sim * (1-rho_sim) * pop15to49) / sum((1-rho_sim) * pop15to49))
  N <- E_recent / (incid * (1-prev) * OmegaT)

  dat <- data
  
  nsamp <- c(rmultinom(1, N, mwsim$pop15to49))
  npos <- rbinom(length(nsamp), nsamp, mwsim$rho_sim)
  nrecent <- rbinom(length(nsamp), npos, mwsim$p_r_i)

  dat$P_i <- npos
  dat$R_i <- nrecent

  return(dat)
}


fit20 <- sampling(hivincid:::stanmodels$incidence, data=dat20[[1]], control = list(adapt_delta = 0.95))

fit200 <- sampling(incidence:::stanmodels$incidence, data=dat200[[2]], control = list(adapt_delta = 0.95))


fit2000 <- sampling(incidence:::stanmodels$incidence, data=dat2000[[2]], control = list(adapt_delta = 0.95))

summary(fit2000, "u_i")$summary[-11*1:2,]
summary(fit200, "u_i")$summary
summary(fit20, "u_i")$summary

#' ## Run simulation on DIDE cluster
workdir <- "/Volumes/jeff/hivincid"
dir.create(workdir)

mrc_config <- didehpc::didehpc_config(workdir=workdir, use_common_lib=FALSE, cluster="mrc", use_workers=TRUE)
pkgsrc <- provisionr::package_sources(local="~/Documents/Research/spatial-estimates/hivincid")
ctx <- context::context_save(file.path(workdir, "context"), packages="hivincid", package_sources=pkgsrc)
mrcq <- didehpc::queue_didehpc(ctx, config=mrc_config)

#' # Generate simulated data
#'
#' 

set.seed(67417478)
Erecent <- c(20, 50, seq(100, 2000, 100))
dat <- lapply(Erecent, function(e) replicate(500, sim_hotspot(e), FALSE))

dat2 <- lapply(dat, lapply, function(d){
  d$Nkappa <- 2;
  d$Xkappa <- cbind(d$Xkappa, as.integer(1:d$N_reg %in% hotid)); d})


tmp <- sampling_summary(hivincid::stanmodels$incidence, data=dat2[[1]][[1]], control=list(adapt_delta=0.95))


mrcq$submit_workers(20)

#' Basic (intercept only) simulation
for(ii in seq_along(Erecent))
  mrcq$mapply(sampling_summary, data=dat[[ii]], name=paste0("bas", Erecent[ii]),
              MoreArgs=list(object=hivincid:::stanmodels$incidence, refresh=0))

for(ii in seq_along(Erecent))
  mrcq$mapply(sampling_summary, data=dat2[[ii]], name=paste0("cov", Erecent[ii]),
              MoreArgs=list(object=hivincid:::stanmodels$incidence, refresh=0,
                            control=list(adapt_delta=0.95)))

mrcq$mapply(sampling_summary, data=dat50, name="sim50", MoreArgs=list(object=hivincid:::stanmodels$incidence, refresh=0))
mrcq$mapply(sampling_summary, data=dat100, name="sim100", MoreArgs=list(object=hivincid:::stanmodels$incidence, refresh=0))
mrcq$mapply(sampling_summary, data=dat200, name="sim200", MoreArgs=list(object=hivincid:::stanmodels$incidence, refresh=0))
mrcq$mapply(sampling_summary, data=dat500, name="sim500", MoreArgs=list(object=hivincid:::stanmodels$incidence, refresh=0))
mrcq$mapply(sampling_summary, data=dat1000, name="sim1000", MoreArgs=list(object=hivincid:::stanmodels$incidence, refresh=0))
mrcq$mapply(sampling_summary, data=dat2000, name="sim2000", MoreArgs=list(object=hivincid:::stanmodels$incidence, refresh=0))
            

#' Collect simulation results
sim20 <- mrcq$task_bundle_get("sim20")$results()
sim50 <- mrcq$task_bundle_get("sim50")$results()
sim100 <- mrcq$task_bundle_get("sim100")$results()
sim200 <- mrcq$task_bundle_get("sim200")$results()
sim500 <- mrcq$task_bundle_get("sim500")$results()
sim1000 <- mrcq$task_bundle_get("sim1000")$results()
sim2000 <- mrcq$task_bundle_get("sim2000")$results()

basall <- list()
for(ii in Erecent){
  print(ii)
  basall[[paste0("bas", ii)]] <- mrcq$task_bundle_get(paste0("bas", ii))$results()
}


covall <- list()
for(ii in Erecent[13:22]){
  print(ii)
  covall[[paste0("cov", ii)]] <- mrcq$task_bundle_get(paste0("cov", ii))$results()
}
                                                      

hotu <- paste0("excess_i[", hotid, "]")
coldu <- paste0("excess_i[", coldid, "]")
allu <- paste0("excess_i[", 1:28, "]")


#' Check Rhat and n_eff

sapply(basall, function(s) mean(apply(sapply(s, "[", allu, "Rhat"), 2, max) < 1.1))
sapply(basall, function(s) mean(apply(sapply(s, "[", allu, "n_eff"), 2, min) < 200))

sapply(covall, function(s) mean(apply(sapply(s, "[", allu, "Rhat"), 2, max) < 1.1))
sapply(covall, function(s) mean(apply(sapply(s, "[", allu, "n_eff"), 2, min) < 200))

#' Remove basulations with Rhat < 1.1

basconv <- lapply(basall, function(s) s[apply(sapply(s, "[", allu, "Rhat"), 2, max) < 1.1])
sapply(basconv, length)
sapply(basconv, function(s) mean(apply(sapply(s, "[", allu, "n_eff"), 2, min) < 200))

covconv <- lapply(covall, function(s) s[apply(sapply(s, "[", allu, "Rhat"), 2, max) < 1.1])
sapply(covconv, length)
sapply(covconv, function(s) mean(apply(sapply(s, "[", allu, "n_eff"), 2, min) < 200))


#' # Probability of identifying hotspots


#' Probability of identifying any hotspots
sapply(basall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, any)))
sapply(covall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, any)))

#' Probability of identifying all hotspots
sapply(basall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, all)))
sapply(covall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, all)))

#' Probability of identifying all hotspots and no coldspots
sapply(basall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, all) &
                                apply(sapply(s, "[", coldu, "2.5%") < 0, 2, all)))
sapply(covall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, all) &
                                apply(sapply(s, "[", coldu, "2.5%") < 0, 2, all)))


hotspots_prob <- data.frame(Erecent = Erecent,
                            bas_any = sapply(basall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, any))),
                            cov_any = sapply(covall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, any))),
                            bas_all = sapply(basall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, all))),
                            cov_all = sapply(covall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, all))),
                            bas_exact = sapply(basall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, all) &
                                                                        apply(sapply(s, "[", coldu, "2.5%") < 0, 2, all))),
                            cov_exact = sapply(covall, function(s) mean(apply(sapply(s, "[", hotu, "2.5%") > 0, 2, all) &
                                                                        apply(sapply(s, "[", coldu, "2.5%") < 0, 2, all))))

devtools::use_data(hotspots_prob)

bas_excess <- sapply(lapply(basall, sapply, "[", allu, "mean"), rowMeans)
cov_excess <- sapply(lapply(covall, sapply, "[", allu, "mean"), rowMeans)

bas_excess <- data.frame(district=data$district, bas_excess)
cov_excess <- data.frame(district=data$district, cov_excess)

devtools::use_data(bas_excess)
devtools::use_data(cov_excess)


#' # PLOTS

devtools::load_all()
library(rstan)
library(maptools)
library(tmap)

data(mwsim)
mwsim@data <- as.data.frame(mwsim@data)


quartz(h=7, w=3.5, pointsize=14)

tm_shape(mwsim) +
  tm_polygons("u_i", palette = "-RdBu", title = "Excess\ntransm.\n(log-scale)",
              style = "cont", colorNA = "white", breaks=seq(-0.75, 0.75, by=0.25), legend.format=list(digits=2, fun=I)) +
  tm_legend(position = c("right","top"))

#' Excess transmission example

set.seed(34164175)
dat20 <- sim_hotspot(20)
dat200 <- sim_hotspot(200)
dat800 <- sim_hotspot(800)

fit20 <- sampling(hivincid:::stanmodels$incidence, data=dat20, control = list(adapt_delta = 0.95))
fit200 <- sampling(hivincid:::stanmodels$incidence, data=dat200, control = list(adapt_delta = 0.95))
fit800 <- sampling(hivincid:::stanmodels$incidence, data=dat800, control = list(adapt_delta = 0.95))

summary(fit20, "excess_i")$summary
summary(fit200, "excess_i")$summary
summary(fit800, "excess_i")$summary


quartz(h=5, w=5.5, pointsize=18)

excess20 <- merge(mwsim, data.frame(district=dat20$district, summary(fit20, "excess_i")$summary[,c("mean", "2.5%", "97.5%")]))

ggplot(data=excess20@data, aes(x=reorder(district, -coordinates(excess20)[,2]))) +
  ggtitle("Excess transmission (log scale)") + 
  geom_abline(intercept = 0, slope=0, colour="grey20", linetype="dashed") +
  geom_point(aes(y=mean), colour="red") +
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=.01, col="red") +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x="District", y="") +
  scale_y_continuous(limits = c(-1.5, 1.5), oob=squish)

quartz(h=6, w=3, pointsize=14)

tm_shape(bas_excess) +
  tm_polygons("bas20",
              palette = "-RdBu", title = "Excess\ntransm.",
              style = "cont", colorNA = "white", breaks=seq(-0.75, 0.75, by=0.25), legend.format=list(digits=2, fun=I)) + 
  tm_layout(panel.labels="E(recent) = 20", panel.label.size = 1.0, legend.outside.size=0.1) +
  tm_legend(position = c("right","top"), title.size=1.0, text.size=0.6, width=0.25, legend.just="right")



#' # E(recent) = 200

quartz(h=5, w=5.5, pointsize=18)

excess200 <- merge(mwsim, data.frame(district=dat200$district, summary(fit200, "excess_i")$summary[,c("mean", "2.5%", "97.5%")]))

ggplot(data=excess200@data, aes(x=reorder(district, -coordinates(excess200)[,2]))) +
  ggtitle("Excess transmission (log scale)") + 
  geom_abline(intercept = 0, slope=0, colour="grey20", linetype="dashed") +
  geom_point(aes(y=mean), colour="red") +
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=.01, col="red") +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x="District", y="") +
  scale_y_continuous(limits = c(-1.5, 1.5), oob=squish)

quartz(h=6, w=3, pointsize=14)

tm_shape(bas_excess) +
  tm_polygons("bas200",
              palette = "-RdBu", title = "Excess\ntransm.",
              style = "cont", colorNA = "white", breaks=seq(-0.75, 0.75, by=0.25), legend.format=list(digits=2, fun=I)) + 
  tm_layout(panel.labels="E(recent) = 200", panel.label.size = 1.0, legend.outside.size=0.1) +
  tm_legend(position = c("right","top"), title.size=1.0, text.size=0.6, width=0.25, legend.just="right")


#' # E(recent) = 800

quartz(h=5, w=5.5, pointsize=18)

excess800 <- merge(mwsim, data.frame(district=dat800$district, summary(fit800, "excess_i")$summary[,c("mean", "2.5%", "97.5%")]))

ggplot(data=excess800@data, aes(x=reorder(district, -coordinates(excess800)[,2]))) +
  ggtitle("Excess transmission (log scale)") + 
  geom_abline(intercept = 0, slope=0, colour="grey20", linetype="dashed") +
  geom_point(aes(y=mean), colour="red") +
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=.01, col="red") +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x="District", y="") +
  scale_y_continuous(limits = c(-1.5, 1.5), oob=squish)

quartz(h=6, w=3, pointsize=14)

tm_shape(bas_excess) +
  tm_polygons("bas800",
              palette = "-RdBu", title = "Excess\ntransm.",
              style = "cont", colorNA = "white", breaks=seq(-0.75, 0.75, by=0.25), legend.format=list(digits=2, fun=I)) + 
  tm_layout(panel.labels="E(recent) = 800", panel.label.size = 1.0, legend.outside.size=0.1) +
  tm_legend(position = c("right","top"), title.size=1.0, text.size=0.6, width=0.25, legend.just="right")



data(bas_excess)

bas_excess <- merge(mwsim, bas_excess)

bas_excess@data[hotid, ]

quartz(h=6, w=12, pointsize=14)

tm_shape(bas_excess) +
  tm_polygons(c("bas20", "bas100", "bas500", "bas1500"),
              palette = "-RdBu", title = "Excess\ntransmission\n(log-scale)",
              style = "cont", colorNA = "white", breaks=seq(-0.75, 0.75, by=0.25), legend.format=list(digits=2, fun=I)) +
  tm_facets(free.scales=FALSE) +
  tm_layout(panel.labels=paste("E(recent) =", c(20, 100, 500, 1500)), panel.label.size = 1.0, legend.outside.size=0.1)
            
  

hotspots_prob_long <- reshape2::melt(hotspots_prob, id="Erecent")

quartz(h=5.5, w=11)

ggplot(subset(hotspots_prob_long, grepl("bas", variable) & Erecent < 1750),
       aes(x=Erecent, y=value, group=variable, color=variable)) +
  geom_point() + geom_line() +
  scale_color_discrete(name="Probability of detecting:",
                         breaks=c("bas_any", "bas_all", "bas_exact"),
                       labels=c("Either hotspot", "Both hotspots", "Only hotspots")) +
  labs(x="E(recent)", y="Probability") +
  theme(text = element_text(size=16))


quartz(h=6, w=7)

ggplot(subset(hotspots_prob_long, Erecent < 1750 & grepl("exact", variable)),
       aes(x=Erecent, y=value, group=variable, color=variable)) +
  geom_point() + geom_line() +
  scale_color_discrete(name="Probability of detecting\ncorrect hotspots:",
                       breaks=c("bas_exact", "cov_exact"),
                       labels=c("Recency\n(& Prev./ART cov)", "Recency +\nCovariate")) +
  labs(x="E(recent)", y="Probability") +
  theme(text = element_text(size=16),
        legend.position="bottom")
