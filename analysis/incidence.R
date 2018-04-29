devtools::load_all(recompile=TRUE)

library(rstan)
library(maptools)
library(data.table)

library(tmap)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

data(mwshest)

mw <- mwshest
mw@data <- data.table(mw@data)

mw$infections <- with(mw@data, pop15to49 * (1-prev_mod5) * prev_mod5 * (1-0.7*artcov_mod5))
mw$pos_dist <- with(mw@data, pop15to49 * prev_mod5)
mw$recent_dist <- with(mw@data, infections / sum(infections))


set.seed(7046752)
mw$nsamp <- c(rmultinom(1, 15255, mw$pop15to49))
mw$npos <- c(rmultinom(1, 1527, mw$pos_dist))
mw$nrecent <- c(rmultinom(1, 18, mw$recent_dist))


mw@data[,list(region, district, recent_dist, nrecent)]
mw@data[,list(sum(recent_dist), nrecent=sum(nrecent), nneg=sum(nsamp - npos)), region]
with(mw@data[,list(sum(recent_dist), nrecent=sum(nrecent), nneg=sum(nsamp - npos)), region], nrecent / nneg * 365/130)



data <- with(mw@data, list(N_reg = nrow(mw),
                           district = mw$district,
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

fit <- sampling(hivincid:::stanmodels$incidence_nat, data=data, control = list(adapt_delta = 0.95))


#' Map showing number of recent cases by district
mw$nrecent_na <- replace(mw$nrecent, mw$nrecent == 0, NA)

quartz(h=7, w=3.5)
tm_shape(mw) +
  tm_polygons("nrecent_na", palette="Reds", title="Number recent", style= "fixed", colorNA="white", breaks=1:4, legend.show = FALSE) +
  tm_text("nrecent", fontface="bold", size=1.2)



print(fit, c("lambda_i", "infections_i", "u_i", "u_raw", "l_rho_i", "l_alpha_i"), include=FALSE, digits_summary=4)


fit <- stan('incidence_nodata.stan', data=data, control = list(adapt_delta = 0.95))
test <- stan('incidence_nodata_test.stan', data=data, control = list(adapt_delta = 0.95))



pars_i <- c(grep("\\_i$", fit5@model_pars, value=TRUE), "u_raw")
print(fit6, pars_i, include=FALSE, digits_summary=4)

with(extract(fit3), plot(colMeans(logit(alpha_anc_i)), colMeans(logit(alpha_i)), pch=19))
with(extract(fit4), plot(colMeans(logit(alpha_anc_i)), colMeans(logit(alpha_i)), pch=19))
with(extract(fit5), plot(colMeans(logit(alpha_anc_i)), colMeans(logit(alpha_i)), pch=19))

with(extract(fit5), plot(colMeans(logit(rho_anc_i)), colMeans(logit(rho_i)), pch=19))

with(mw, points(logit(ancrt_pos / ancrt_n), logit(npos / nsamp), col=2, pch=19))


fit0 <- stan('incidence.stan', data=data, control = list(adapt_delta = 0.95))

print(fit0, c("lambda_i", "infections_i", "u_i", "u_raw"), include=FALSE, digits_summary=4)

post <- extract(fit)

mw$u_i <- colMeans(post$u_i)

save(fit, data, mw, file="~/Documents/Meetings/2018/Jan - Mar/01. Imperial Statistics seminar/incidence-outputs.RData")

## Create adjacency matrix
adj <- spdep::poly2nb(mw)
## adj <- adj[-28]  # Remove Likoma
n_nb <- spdep::card(adj)

adj <- data.frame(reg_i=rep(seq_along(adj), n_nb+1), reg_j=unlist(mapply(c, seq_along(adj), adj)))

do.call(rbind, mapply(data.frame, i=seq_along(adj), j=adj, SIMPLIFY=FALSE))


#' # Plot results

quartz(h=7, w=3.5, pointsize=14)


mwfit2 <- merge(mw, dcast(estall[model == "Case 2" & outcome == "mean"], district ~ param))


mwfit2 <- merge(mwfit, setNames(data.frame(data$district, summary(fit, "rho_i")$summary[,c("mean", "2.5%", "97.5%")]), c("district", "rho", "rho_cil", "rho_ciu")))
mwfit <- merge(mwfit, setNames(data.frame(data$district, summary(fit, "alpha_i")$summary[,c("mean", "2.5%", "97.5%")]), c("district", "alpha", "alpha_cil", "alpha_ciu")))
mwfit <- merge(mwfit, setNames(data.frame(data$district, summary(fit, "u_i")$summary[,c("mean", "2.5%", "97.5%")]), c("district", "u_i", "u_i_cil", "u_i_ciu")))
mwfit <- merge(mwfit, setNames(data.frame(data$district, summary(fit, "lambda_i")$summary[,c("mean", "2.5%", "97.5%")]), c("district", "lambda", "lambda_cil", "lambda_ciu")))


## Map of prevalence, ART, coverage, u_i, and incidence rate by district
quartz(h=6, w=2.5, pointsize=14)

tm_shape(mwfit) +
  tm_polygons("u_i", palette = "Greens", style = "cont", colorNA = "white",
              breaks=seq(-0.75, 0.75, by=0.25), legend.format=list(digits=2, fun=I),
              title="") +
  tm_legend(position = c("right","top")) + 
  tm_layout(panel.labels=expression(u[i]))

tm_shape(mwfit) +
  tm_polygons("rho", palette = "Purples", style = "cont", colorNA = "white",
              breaks=seq(0.0, 0.3, 0.05),
              legend.format=list(digits=2, fun=I),
              title="") +
  tm_legend(position = c("right","top")) + 
  tm_layout(panel.labels="Prevalence")

tm_shape(mwfit) +
  tm_polygons("alpha", palette = "Blues", style = "cont", colorNA = "white",
              breaks=seq(0.4, 0.8, 0.1),
              legend.format=list(digits=2, fun=I),
              title="") +
  tm_legend(position = c("right","top")) + 
  tm_layout(panel.labels="ART coverage")

mwfit$loglambda <- log(mwfit$lambda)

tm_shape(mwfit) +
  tm_polygons("loglambda", palette = "-Reds", style = "cont", colorNA = "white",
              breaks=log(seq(0.001, 0.011, 0.002)),
              legend.format=list(fun=exp),
              title="") +
  tm_legend(position = c("right","top")) + 
  tm_layout(panel.labels="Incidence rate", legend.format=list(digits=3))




## Vertical bars of incidence rate by district
quartz(h=4, w=8, pointsize=18)

ggplot(data=mwfit@data, aes(x=reorder(district, -coordinates(mwfit)[,2]))) +
  ggtitle("District-level incidence rate") + 
  geom_abline(intercept = summary(fit, "lambda")$summary[1], slope=0, colour="grey20", linetype="dashed") +
  geom_point(aes(y=lambda), colour="red") +
  geom_errorbar(aes(ymin=lambda_cil, ymax=lambda_ciu), width=.01, col="red") +
  theme(axis.text.x = element_text(angle=60, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  labs(x="District", y="Incidence rate") +
  scale_y_continuous(limits = c(0, 0.015), oob=squish)
