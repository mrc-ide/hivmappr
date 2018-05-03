#' # Introduction
#'
#' This script outlines the steps for a shiny app workflow with example data. Much of this needs to be made into functions.
#'

devtools::load_all()
library(rgdal)
library(magrittr)
library(rstan)
library(data.table)
library(RColorBrewer)

library(ggplot2)
library(ggridges)
library(gridExtra)


#' ## 1. Upload shape file, generate a spreadsheet table to input data

## Shape file from DHS: https://spatialdata.dhsprogram.com/boundaries/#view=map&countryId=MW&surveyId=483&level=2
sh <- readOGR(system.file("extdata", "mwsh", package="hivincid"))
plot(sh)

#' ## 2. Copy and paste data into spreadsheet

df <- read.csv(system.file("extdata", "mwdf.csv", package="hivincid"))
mw <- merge(sh, df)


#' ## 3. Do some basic data visualization and checks 

## TO BE COMPLETED

#' ## 4. Press a button to fit the model

#' Prepare data for Stan model input
data <- with(mw@data, list(N_reg = length(district),
                           district = district,
                           prev_est = prev_survey,
                           prev_se = prev_survey_se,
                           anc1_obs = cbind(ancrt_neg = ancrt_n - ancrt_pos,
                                            ancrt_noart=ancrt_pos - ancrt_art,
                                            ancrt_art=ancrt_art),
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
                           Xkappa = matrix(10, length(district), 1),
                           sigma_u_sd = 1))


#' Fit stan model
fit <- sampling(stanmodels$incidence_rita, data=data, control = list(adapt_delta = 0.95))

#' ## 5. Display maps and posterior figure

#' ### Map of posterior mean estimates

## extract summary estimates 

est <- summary(fit, c("rho_i", "alpha_i", "u_i", "lambda_i", "infections_i"))$summary[, "mean"]
est <- data.table(param = names(est), value=est)
est$district_idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", est$param))
est$district <- factor(mw$district[est$district_idx], data$district)
est$region <- factor(mw$region[est$district_idx], c("Northern", "Central", "Southern"))
est$param <- sub("([^\\[]+).*", "\\1", est$param)


## Extract map polygon data 
mwpoly <- map_data(mw, namefield="district")

th <- function(title)
  list(geom_map(aes(fill = value), map = mwpoly),
       expand_limits(x = mwpoly$long, y = mwpoly$lat),
       labs(x = NULL, y = NULL, title=title),
       coord_map(),
       theme(axis.ticks = element_blank(),
             axis.text = element_blank(),
             panel.background = element_blank(),
             panel.border = element_blank(),
             legend.position=c(1, 1),
             legend.justification=c(0.5, 1),
             legend.background = element_rect(fill = NA),
             legend.key.size=unit(12, "points"),
             legend.text = element_text(size=10),
             plot.title=element_text(hjust=0.5, face="bold"),
             text = element_text(size=10)))

panA <- ggplot(est[param == "rho_i"], aes(map_id = district)) +
  scale_fill_distiller(element_blank(), palette="Purples", direction=1,
                       lim=c(0, 0.25), labels=function(x) round(100*x)) +
  th("Prevalence (%)")
panB <- ggplot(est[param == "alpha_i"], aes(map_id = district)) +
  scale_fill_distiller(element_blank(), palette="Blues", direction=1,
                       lim=c(0.3, 0.8), labels=function(x) round(100*x)) +
  th("ART Coverage (%)")
panC <- ggplot(est[param == "u_i"], aes(map_id = district)) +
  scale_fill_distiller(element_blank(), palette="PRGn", direction=1,
                       lim=c(-0.2, 0.2)) +
  th("u_i")
panD <- ggplot(est[param == "lambda_i"], aes(map_id = district)) +
  scale_fill_distiller(element_blank(), palette="Reds", direction=1,
                       lim=c(1, 8)/1e3, labels=function(x) round(1000*x)) +
  th("Incidence / 1000")  

quartz(h=3, w=6)
grid.arrange(panA, panB, panC, panD, ncol=4)


#' ### Density plot of posterior district estimates

samp <- as.matrix(fit, c("rho_i", "alpha_i", "u_i", "lambda_i", "infections_i")) %>%
  abind::abind(along=0)
names(dimnames(samp)) <- c("model", NA, "param")
samp <- melt(samp) %>% data.table
samp$district_idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", samp$param))
samp$district <- factor(mw$district[samp$district_idx],
                        rev(mw$district))
samp$region <- factor(mw$region[samp$district_idx], c("Northern", "Central", "Southern"))
samp$param <- sub("([^\\[]+).*", "\\1", samp$param)

th <- list(theme_minimal(),
           theme(plot.title=element_text(hjust=0.5, face="bold", size=8),
                 text = element_text(size=10),
                 aspect.ratio=2,
                 axis.text.y=element_blank(),
                 axis.ticks.y = element_line(),
                 axis.ticks.length = unit(2, "pt")))

panA <- ggplot(data=samp[param == "rho_i" & model == "rita"],
               aes(x=value, y=district, fill = ..x.. )) +
  geom_density_ridges_gradient(rel_min_height=0.01) +
  geom_vline(xintercept=0.10, color="grey20", linetype="dashed") +
  scale_x_continuous(element_blank(), labels=function(x) round(100*x, 1)) +
  scale_y_discrete(element_blank()) +
  scale_fill_distiller(guide = "none", palette="Purples", direction=1, trans="log10") +
  labs(title="Prevalence (%)") +
  th + theme(axis.text.y=element_text(size=7, hjust=1))

panB <- ggplot(data=samp[param == "alpha_i" & model == "rita"],
               aes(x=value, y=district, fill = ..x.. )) +
  geom_density_ridges_gradient(rel_min_height=0.01) +
  geom_vline(xintercept=0.51, color="grey20", linetype="dashed") +
  scale_x_continuous(element_blank(), labels=function(x) round(100*x, 1)) +
  scale_y_discrete(element_blank()) +
  scale_fill_distiller(guide = "none", palette="Blues", direction=1) +
  labs(title="ART Coverage (%)") +
  th

panC <- ggplot(data=samp[param == "lambda_i" & model == "rita"],
       aes(x=value, y=district, fill = ..x.. )) +
  geom_density_ridges_gradient(rel_min_height=0.025) +
  geom_vline(xintercept=0.0036, color="grey20", linetype="dashed") +
  scale_x_log10(element_blank(),
                breaks=c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02),
                limits=c(0.0006, 0.025), labels=function(x) round(1000*x, 1)) +
  scale_y_discrete(element_blank()) +
  scale_fill_distiller(guide = "none", palette="Reds", direction=1) +
  labs(title="Incidence / 1000 (log)") +
  th

panD <- ggplot(data=samp[param == "infections_i" & model == "rita"],
       aes(x=value, y=district, fill = ..x.. )) +
  geom_density_ridges_gradient(rel_min_height=0.025) +
  scale_x_log10(element_blank(), limits=c(50, 10000), labels=function(x) round(x)) +
  scale_y_discrete(element_blank()) +
  scale_fill_distiller(guide = "none", palette="Reds", direction=1) +
  labs(title="New infections", fontface="bold") +
  th

quartz(h=3, w=6)
grid.arrange(panA, panB, panC, panD, ncol=4, widths=c(1.35, 1, 1, 1))

#' ### Something like the violin plots showing data inputs vs. estimates

## TO BE COMPLETED

#' ## 6. Produce table of outputs to view and download

est <- summary(fit, c("rho_i", "alpha_i", "lambda_i", "infections_i"))$summary[, c("mean", "2.5%", "97.5%")]
colnames(est) <- c("mean", "ci_l", "ci_u")
est <- data.table(param = rownames(est), est)
est$district_idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", est$param))
est$district <- factor(mw$district[est$district_idx], data$district)
est$region <- factor(mw$region[est$district_idx], c("Northern", "Central", "Southern"))
est$param <- sub("([^\\[]+).*", "\\1", est$param) %>% sub("\\_i$", "", .)
est$scale <- c(rho = 100, alpha = 100, lambda = 1000, infections = 1)[est$param]
est$digits <- c(rho = 1, alpha = 0, lambda = 1, infections = -2)[est$param]
est$str <- with(est, sprintf(paste0("%.", pmax(digits, 0), "f (%.", pmax(digits, 0), "f, %.", pmax(digits, 0), "f)"),
                             round(scale*mean, digits),
                             round(scale*ci_l, digits),
                             round(scale*ci_u, digits)))
est$label <- factor(est$param, c("rho", "alpha", "lambda", "infections"),
                    c("Prevalence (%)", "ART coverage (%)", "Incidence (per 1000)", "New infections"))

dcast(est, district+region ~ label, value.var="str") %>% knitr::kable()
