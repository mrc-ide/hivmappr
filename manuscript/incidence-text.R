#' ---
#' title: A small-area model for estimating subnational HIV incidence
#' author: Jeffrey W Eaton, Sumali Bajaj, \textit{and others}
#' output: 
#'  pdf_document:
#'    number_sections: true
#' ---

##+ setup, include=FALSE
library(knitr)
opts_chunk$set(tidy=TRUE, warning=FALSE, cache=TRUE, message=FALSE)
devtools::load_all()
library(magrittr)
library(rstan)
library(data.table)
library(RColorBrewer)

library(ggplot2)
library(ggridges)
library(gridExtra)

#' 
#' # Introduction
#' 
#' We have previously described a model for estimating subnational HIV prevalence and ART coverage
#' using survey data and routine health system data [@Eaton2017]. This note describes an
#' extension of our model to estimate subnational HIV incidence including probabilistic
#' uncertainty. We illustrate the proposed model for estimating the number of new HIV infections
#' by subnational area for two-common cases of data availability:
#'
#' 1. Subnational data about HIV prevalence and ART coverage are available, but estimates of HIV
#' incidence are only available at the national level, for example from national EPP/Spectrum
#' estimates.
#' 1. Direct data about subnational HIV incidence are available through inclusion of tests
#' for recent infection in population survey data.
#'
#' Finally, we discuss the model as a framework for identifying HIV transmission `hotspots' for
#' prioritising further HIV prevention activities.
#'
#' 
#' # Modelling HIV incidence
#'
#' The key principle underpinning the modelling of infectious disease dynamics is that the risk of
#' acquiring infection depends on (1) the rate of tranmsission from an infectious host and (2)
#' the probability of contacting an infectious host [@Anderson1982]. For HIV, it is useful to express
#' the HIV incidence rate, or *force of infection*, $\lambda$ as the product of the transmission rate
#' $\kappa$ and the prevalence of unsuppressed HIV viral load in the population, that is
#' $$\lambda = \kappa \cdot \rho \cdot (1-\omega\cdot\alpha),$$
#' where $\rho$ is the HIV prevalence, $\alpha$ is the ART coverage, and $\omega$ is the relative
#' reduction in tranmission for persons on ART. This relationship between HIV incidence and prevalence
#' has long been used to estimate and project national HIV incidence trends with the EPP model [@EPP].
#'
#' This theory suggests a log-linear model for the HIV incidence rate $\lambda_i$ in region $i$
#' \begin{equation}
#'   \log(\lambda_i) = \log(\kappa_0) + \log(\rho_i) + \log(1-\omega\cdot\alpha_i) + u_i,
#' \end{equation}
#' where $\kappa_0$ is the average transmission rate, $\rho_i$ is the HIV prevalence
#' in region $i$, $\alpha_i$ is the ART coverage, and $u_i$ is a random effect capturing other
#' region-level differences in the relative transmission rate. Removing the prevalence term from the linear
#' model to directly model spatial patterns in the transmission rate $\kappa_i$ may be more intuitive:
#' \begin{equation}
#'   \log(\kappa_i) = \log(\kappa_0) + \log(1-\omega\cdot\alpha_i) + u_i,
#' \end{equation}
#'
#' In this formulation, the incidence rate $\lambda_i = \kappa_i\cdot\rho_i$, noting that
#' $\kappa_i$ is the ratio of incidence over prevalence. The expression for the annual number of new
#' infections $I_i$ in region $i$ is
#' $$I_i = \kappa_i \cdot \rho_i \cdot ((1-\rho_i) \cdot N_i),$$
#' where $N_i$ is the total population and observing that $(1-\rho_i) \cdot N_i$ is the number
#' susceptible in the population.
#'
#' ## Case 1: No directly observed HIV incidence
#'
#' In many cases, data about HIV prevalence and ART coverage are available at subnational levels from population
#' surveys and routine programme data, but estimates for HIV incidence are only available at national or first
#' administrative unit level, for example from the EPP model. In this case, the objective for a model for
#' subnational HIV incidence is to estimate the most plausible subnational distribution of new HIV cases based on
#' the available information about subnational prevalence and ART coverage $\rho_i$ and $\alpha_i$, and such that
#' the aggregated subnational estimates of new infections are consistent with the exogenously estimated national
#' (or larger area) esimates. The random effect terms $u_i$ allow for additional uncertainty about subnational new
#' infections arising from transmission being higher or lower than predicted solely by the prevalence of unsuppressed
#' viral load.
#'
#' Exogenous estimates and uncertainty about the national HIV incidence ratey $\lambda^{Nat}$ may be summarized as a
#' normal distribution for $\log(\lambda^{Nat})$:
#' \begin{equation}
#'   \label{eqn:incid-prior}
#'   \log(\lambda^{Nat})\sim \mathrm{Normal}(\log(\hat\lambda^{Nat}), \sigma_{\lambda^{Nat}})
#' \end{equation}
#' 
#' We now seek to relate the relate the national HIV incidence rate to the aggregation of the subnational incidnece rates.
#' Defining $S_i=(1-\rho_i)\cdot N_i$ as the number not infected with (susceptible to) HIV infection in region $i$,
#' the national HIV incidence rate is expressed as
#' \begin{equation}
#'   \lambda^{Nat} = \frac{\sum_i \lambda_i \cdot S_i}{\sum_i S_i}
#'       = \frac{\sum_i \kappa_0 \cdot(1-\omega \cdot \alpha_i) \cdot \rho_i \cdot \exp(u_i) \cdot S_i}{\sum_i S_i}.
#' \end{equation}
#'
#' Taking the logarithm and rearranging terms yields that
#' \begin{equation}
#' \label{eqn:loglambdanat}
#'   \log(\lambda^{Nat}) = \log(\kappa_0) + \log\left(\sum_i (1-\omega \cdot \alpha_i) \cdot \rho_i \cdot \exp(u_i) \cdot S_i\right) - \log\left(\sum_i S_i\right)
#' \end{equation}
#'
#' Substituting equation (\ref{eqn:loglambdanat}) into equation (\ref{eqn:incid-prior}) defines a model for subnational incidence given estimates of prevalence, ART coverage, and the variance of the spatial random effects:
#' \begin{equation}
#'   \begin{aligned}
#'     u_i &\sim \mathrm{Normal}(0, \sigma_u) \\
#'     \log(\kappa_0) | \rho_i, \alpha_i, u_i &\sim \mathrm{Normal}\left(\log(\hat\lambda^{Nat}) - \log\left(\sum_i k_i \cdot S_i\right) + \log\left(\sum_i S_i\right), \sigma_{\lambda^{Nat}}\right) \\
#'     \lambda_i &= \kappa_0 \cdot k_i,
#'   \end{aligned}
#' \end{equation}
#' where $k_i = \rho_i \cdot (1-\omega \cdot \alpha_i) \cdot \exp(u_i)$. The joint distribution for ${\rho_i, \alpha_i}$ is estimated using the model we have described previously for subnational HIV prevalence and ART coverage. In the absence of subnational HIV incidence data, a value for $\sigma_u$ must be assumed or derived from other sources.
#' 
#' 
#' ## Case 2: Recent infection testing algorithm
#'
#' A number of national HIV household surveys are now including biomarker-based laboratory algorithms
#' to identify whether HIV infections were `recently' acquired. Instead of simply providing estimates
#' of the number HIV positve and number HIV negative in the population, these surveys furnish a
#' trinomial observation $\{R_i, L_i, N_i\}$ of the number HIV positive and recently infected, the
#' number of long-term (not-recent) HIV infected, and the number HIV negative, respectively, in each
#' region $i$.
#'
#' Kassanjee and colleagues describe an estimator for HIV incidence as a function of
#' the HIV prevalence $\rho$, the proportion of HIV-positive cases that were recently infected $p^R$,
#' the mean duration of recent infection (MDRI) $\Omega_T$, the expected proportion of long-term
#' infections misclassified as recent $\beta_T$, and the recency cut-off period $T$:
#' \begin{equation}
#'   \label{eq:kassanjee}
#'   \lambda = \frac{(p^R - \beta_T)\cdot \rho}{(1-\rho)\cdot(\Omega_T - \beta_TT)}.
#' \end{equation}
#'
#' Rearranging equation \ref{eq:kassanjee}, the expected proportion of recent infections $p^R_i$ in a
#' region as a function of the prevalence $\rho_i$ and force of infection $\lambda_i$ is
#' \begin{equation}
#' \begin{aligned}
#' p^R_i &= \frac{\lambda_i \cdot (1-\rho_i)\cdot(\Omega_T - \beta_TT) + \beta_T\rho_i}{\rho_i} \\
#'  &= \kappa_i \cdot (1-\rho_i)\cdot(\Omega_T - \beta_TT) + \beta_T
#' \end{aligned}
#' \end{equation}
#' where the second equality follows from substituting $\lambda_i = \kappa_i\cdot\rho_i$.
#'
#' Using this expression, we may express a likelihood for the number recently $R_i$
#' of the total number HIV positive $P_i = R_i + L_i$:
#' \begin{equation}
#'   \label{eqn:recent-lik}
#'   R_i \sim \mathrm{Binomial}(p^R_i, P_i)
#' \end{equation}
#'
#' By convention, we assume
#' \begin{equation}
#'   \begin{aligned}
#'     \Omega_T &\sim \mathrm{Normal}(\Omega_{T0}, \sigma^2_{\Omega}) \\
#'     \beta_T &\sim \mathrm{Normal}_{\beta_T\geq0}(\beta_{T0}, \sigma^2_\beta),
#'   \end{aligned}
#' \end{equation}
#' where the values of $\Omega_{T0}$, $\sigma^2_{\Omega}$, $\beta_{T0}$, $\sigma^2_\beta$, and $T$
#' are specified for a given survey depending on the HIV subtypes in the survey population and
#' the details of the RITA used.
#'
#' In realistic applications of the model, data $R_i$ and $P_i$ most likely arise from complex survey
#' data including clustered sampling with unequal weights, which must be accounted for in analysis. We
#' use normalized weighted counts for $R_i$ and $P_i$ for evaluation of the likelihood in equation~
#' (\ref{eqn:recent-lik}). A conventional approach for model-based analysis of survey data is to
#' approximate the likelihood for the observed counts with design-based estimates and standard errors
#' at the level of analytical interest, in this case a direct estimate of $\hat p^R_i$. This approach is
#' infeasible here because recent infection is a rare event, and indeed we expect a large number of
#' regions with $R_i = 0$.  The use of weighted counts ensures that aggregated model-based
#' estimates are consistent with design-based estimates, though may understate the uncertainty arising
#' from the clustered sampling \footnote{This approach is consistent with the primary analysis of PHIA
#' survey results which assumed a design effect of 1.0 for the proportion recent in published incidence
#' estimates.}. The bias will be minimal if the regions $i$ are small enough that sampled clusters within
#' a given region are relatively homogeneous. Optimal approaches for analysis of complex survey data with
#' small counts remains an open research question. Modelling the full survey design, for example by
#' post-stratifying the analysis by all factors used in defining analysis weights and including
#' cluster random effects nested with area-level effects, has been used elsewhere, but is overly
#' cumbersome for our exposition and may be impossible if full survey design details are not avaialble to
#' the analyst.
#' 
#' # Application to simulated data from Malawi
#'
#' We demonstrate the models for subnational incidence estimates using _simulated data_ representative of
#' Malawi. National-level estimates for HIV prevalence, ART coverage, and incidence are available from the
#' Malawi Population HIV Impact Assessment (MPHIA) survey conducted in 2016. We generate simulated subnational
#' data consistent with national-level survey results for illustration. Further details of the simulated data
#' are described in the Appendix and the simulated district-level survey results are presented in
#' Table~\ref{tab:}.
#'
#' For `case 1' in which only national incidence estimates are available, we used estimates for national HIV
#' incidence among adults 15-49 years of 0.36\% (SE 0.087\%). Log-transforming these estimates yields estimates
#' $\log(\hat\lambda^{Nat}) = -5.63$ and $\sigma_{\lambda^{nat}} = 0.241$. These 
#'
#' Throughout all analyses, we used the value $\omega=0.7$ as assumed by the EPP model for the average
#' reduction in transmission per percentage increase in ART coverage
#'

##+ simulate data, include=FALSE
set.seed(7046752)
data(mwshest)
mw <- mwshest
mw@data[c("lat", "long")] <- coordinates(mw)

omega <- 0.7

## Number HIV positive in each district
mw$hivpos <- with(mw@data, pop15to49 * prev_mod5)

## Distribution of HIV infections by district, determined by population size,
## prevalence of unsuppressed viral load, and number susceptible.
mw$infections <- with(mw@data, hivpos * (1-prev_mod5) * (1-0.7*artcov_mod5))

## Distribute number sampled in each district proportional to population size
mw$nsamp <- c(rmultinom(1, 15255, mw$pop15to49))

## Sample number HIV+ proportional to distribution of HIV population
mw$npos <- c(rmultinom(1, 1527, mw$hivpos))

## Sample 18 observed recent cases proporitonal to new infections
mw$nrecent <- c(rmultinom(1, 18, mw$infections))

## Survey prevalence estimates and standard error, assume design effect = 2.0
mw$prev_est <- with(mw@data, npos / nsamp)
mw$prev_se <- with(mw@data, sqrt(2 * prev_est * (1 - prev_est) / nsamp))

##+ model fitting, include = FALSE

## create index to order by region then north to south
distord <- order(factor(mw$region, c("Northern", "Central", "Southern")), -coordinates(mw)[,2])

data <- with(mw@data, list(N_reg = nrow(mw),
                           district = mw$district,
                           prev_est = prev_est,
                           prev_se = prev_se,
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
                           Xkappa = matrix(10, length(district), 1),
                           log_lambda_nat_mean = log(0.0036),
                           log_lambda_nat_sd = (0.0053 - 0.0019)/(2*qnorm(0.975)) / 0.0036,
                           sigma_u_sd = 1))


## Model fit for case 1
fit1 <- sampling(stanmodels$incidence_nat, data=data, control = list(adapt_delta = 0.99))

## Model fit for case 2
fit2 <- sampling(stanmodels$incidence_rita, data=data, control = list(adapt_delta = 0.99))

fitall <- list(nat=fit1, rita=fit2)

#' The 
##+ map, echo=FALSE, results="hide", fig.height=3, fig.width=6

## Extract district prevalence and ART coverage estimates for summarizing and mapping

est1 <- summary(fit1, c("rho_i", "alpha_i", "u_i", "lambda_i", "infections_i"))$summary
est2 <- summary(fit2, c("rho_i", "alpha_i", "u_i", "lambda_i", "infections_i"))$summary
names(dimnames(est1)) <- names(dimnames(est2)) <- c("param", "outcome")

estall <- rbind(data.frame(model="nat", melt(est1)),
                data.frame(model="rita", melt(est2))) %>% data.table()
estall$district_idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", estall$param))
estall$district <- factor(mw$district[estall$district_idx], mw$district[distord])
estall$region <- factor(mw$region[estall$district_idx], c("Northern", "Central", "Southern"))
estall$param <- sub("([^\\[]+).*", "\\1", estall$param)
estmean <- estall[model == "rita" & outcome == "mean"]

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

panA <- ggplot(estmean[param == "rho_i"], aes(map_id = district)) +
  scale_fill_distiller(element_blank(), palette="Purples", direction=1,
                       lim=c(0, 0.25), labels=function(x) round(100*x)) +
  th("Prevalence (%)")
panB <- ggplot(estmean[param == "alpha_i"], aes(map_id = district)) +
  scale_fill_distiller(element_blank(), palette="Blues", direction=1,
                       lim=c(0.3, 0.8), labels=function(x) round(100*x)) +
  th("ART Coverage (%)")
panC <- ggplot(estmean[param == "u_i"], aes(map_id = district)) +
  scale_fill_distiller(element_blank(), palette="PRGn", direction=1,
                       lim=c(-0.2, 0.2)) +
  th("u_i") +
  geom_text(data=mw@data, aes(x=lat, y=long, label=nrecent),
            fontface="bold", col="grey25", size=2.5)
panD <- ggplot(estmean[param == "lambda_i"], aes(map_id = district)) +
  scale_fill_distiller(element_blank(), palette="Reds", direction=1,
                       lim=c(1, 8)/1e3, labels=function(x) round(1000*x)) +
  th("Incidence / 1000")  

grid.arrange(panA, panB, panC, panD, ncol=4)

#' The next figure shows the posterior distribution for HIV prevalence, ART coverage,
#' HIV incidence rate, and the number of new infections by district. Vertical
#' dashed lines in the first three panels demarcate the national estimate for each
#' outcome. The figure illustrates the large heterogeneity in HIV prevalence across
#' districts and district-level prevalence estimates are relatively precise, reflecting
#' the large amount of available data about HIV prevalence. Estimates for the HIV incidence
#' rate are more uncertain given that district estimates are based on a total of 18 cases
#' of recent infection, but overall the model estimates higher incidence rate in high
#' prevalence districts in Southern Malawi. The predicted number of new infections
#' reflects the variation in incidence rate and the district population. The largely
#' urban districts of Blantyre and Lilongwe account for an estimated 30\% (21\%--42\%)
#' of all new infections. Blantyre is one of the highest prevalence and incidence
#' districts, but prevalence and incidence in Lilongwe are estimated to be slightly
#' lower than the national levels.
#' 

##+ ggridges rita, echo=FALSE, results="hide", fig.height=3, fig.width=6

sampall <- fitall %>% lapply(as.matrix, c("rho_i", "alpha_i", "u_i", "lambda_i", "infections_i"))

sampall <- abind::abind(sampall, along=0)
names(dimnames(sampall)) <- c("model", NA, "param")
sampall <- melt(sampall)
sampall$district_idx <- as.integer(sub(".*\\[([0-9]+)\\]", "\\1", sampall$param))
sampall$district <- factor(mw$district[sampall$district_idx],
                           rev(mw$district[distord]))
sampall$region <- factor(mw$region[sampall$district_idx], c("Northern", "Central", "Southern"))
sampall$param <- sub("([^\\[]+).*", "\\1", sampall$param)
sampall <- data.table(sampall)

## Proportion of all new infections in Lilongwe and Blantyre 
inf_urb <- rowSums(as.matrix(fit2, c("infections_i[10]", "infections_i[18]")))
inf_all <- as.matrix(fit2, "infections")
mean(inf_urb / inf_all)
quantile(inf_urb / inf_all, c(0.025, 0.975))

th <- list(theme_minimal(),
           theme(plot.title=element_text(hjust=0.5, face="bold", size=8),
                 text = element_text(size=10),
                 aspect.ratio=2,
                 axis.text.y=element_blank(),
                 axis.ticks.y = element_line(),
                 axis.ticks.length = unit(2, "pt")))

panA <- ggplot(data=sampall[param == "rho_i" & model == "rita"],
               aes(x=value, y=district, fill = ..x.. )) +
  geom_density_ridges_gradient(rel_min_height=0.01) +
  geom_vline(xintercept=0.10, color="grey20", linetype="dashed") +
  scale_x_continuous(element_blank(), labels=function(x) round(100*x, 1)) +
  scale_y_discrete(element_blank()) +
  scale_fill_distiller(guide = "none", palette="Purples", direction=1, trans="log10") +
  labs(title="Prevalence (%)") +
  th + theme(axis.text.y=element_text(size=7, hjust=1))
####
####
panB <- ggplot(data=sampall[param == "alpha_i" & model == "rita"],
               aes(x=value, y=district, fill = ..x.. )) +
  geom_density_ridges_gradient(rel_min_height=0.01) +
  geom_vline(xintercept=0.51, color="grey20", linetype="dashed") +
  scale_x_continuous(element_blank(), labels=function(x) round(100*x, 1)) +
  scale_y_discrete(element_blank()) +
  scale_fill_distiller(guide = "none", palette="Blues", direction=1) +
  labs(title="ART Coverage (%)") +
  th
####
####
panC <- ggplot(data=sampall[param == "lambda_i" & model == "rita"],
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
####
####
panD <- ggplot(data=sampall[param == "infections_i" & model == "rita"],
       aes(x=value, y=district, fill = ..x.. )) +
  geom_density_ridges_gradient(rel_min_height=0.025) +
  scale_x_log10(element_blank(), limits=c(50, 10000), labels=function(x) round(x)) +
  scale_y_discrete(element_blank()) +
  scale_fill_distiller(guide = "none", palette="Reds", direction=1) +
  labs(title="New infections", fontface="bold") +
  th

grid.arrange(panA, panB, panC, panD, ncol=4, widths=c(1.35, 1, 1, 1))


#' TODO: Compare posterior estimates for national prior vs. RITA data.
#'
#' TODO: Compare incidence estimates depending on assumption about standard deviation of $u_i$.


#'
#' # Identifying HIV transmission hotspots
#'
#' 
#' # Dicussion and future directions
#'
#' We have proposed a small-area model for subnational HIV incidence that combines basic theory of
#' infectious disease transmission and a stochastic model for the force of infection. We use the
#' prevalence of unsuppressed HIV viral load to predict the relative levels of regional HIV incidence.
#' Including this theoritical epidemiologic relationship in the model is important because in many
#' cases there are not directly observed data about HIV incidence, and where data are available,
#' the number of recent infections observed in subnational regions is small, in many cases zero.
#' We anticipate that validation of this assumption through standard statistical approaches such
#' as cross-validation may be challenging due to the paucity of subnational HIV incidence data.
#' However, the assumption is supported by recent studies from Kwa-Zulu Natal South Africa and
#' western Kenya that found that HIV incidence was highly correlated with the prevalence of
#' unsuppressed HIV viral load in the loca area [@Vandormael2017, @Ndhiwa abstract].
#'
#' The model for $\log(\lambda_i)$ may be enhanced in natural ways, including adding fixed effects
#' that predict regions of higher or lower HIV transmission, or through more sophisticated modelling
#' of the $u_i$, such as capturing spatial correlation. Fixed effects may include factors such as urbanization,
#' distance from major roadways, or locations of risk enumerated through risk-mapping exercises. With
#' sufficient data on subnational HIV incidence, these relationships may be learned through regression,
#' or may be included in the model through informative prior distributions derived from epidemiological
#' literature.
#'
#' One appealing feature of the the proposed small-area incidence model is the interpretation of
#' the $u_i$ random-effect terms. Positive or negative values of $u_i$ indicate regions where HIV
#' transmission is higher or lower than would be expected based on the local prevalence of
#' unsuppressed HIV viral alone. This presents a framework for statistically identifying HIV transmission
#' 'hotspots', which may be useful for focusing HIV prevention programmes for greater impact. Given the
#' small number of recent infection cases expected in any one subnational area, we do not anticipate that
#' analysis of recent infection testing algorithms in household surveys is likely to provide sufficient
#' statistical power to robustly identify transmission hotspots. However, more precise estimates for
#' transmission hotspots may be possible through extending the model to incorporate other data about
#' or covariates for HIV, such as case reports of new HIV diagnoses from routine programmatic data or
#' estimates of HIV incidence derived from self-reported HIV testing history in household surveys
#' [@Fellows].
#'
#' ## Relationship between the force of infection and prevalence
#'
#' - Expect higher force of infection in areas of higher prevalence.
#' - Not necessasrily true.
#' - Including $\beta\cdot\log(\rho_i)$ regression term in model for $\lambda_i$, hypothesise positive coefficient
#'
#' ## Indirect inference about incidence from prevalence trends
#'
#' Noting that the proposed small-area HIV incidence model for uses the exact same functional form for
#' predicting HIV incidence as a function of the force of inection, prevalence, and ART coverage as the
#' EPP model uses for inferring HIV incidence trends, we begin to see that we are conceptually only a
#' small step from a in integrated model for jointly estimating spatio-temporal HIV prevalence, incidence,
#' and service coverage. The final component requires specifying a model for HIV mortality, such that
#' we can model the relationship between HIV incidence, prevalence, and mortality over time.
#' However, acheiving this is anticipated to entail substantial research effort both in model development
#' and computational aspects.
#'
#' \appendix
#' -------------------
#' # Appendix: simulated subnational survey data for Malawi
#'
#' The 2016 MPHIA survey included RITA-based incidence estimates. Regional data from this survey are not
#' available for analysis, but it has been reported that among respondents aged 15-49 at the national level,
#' there were 13727.9 HIV-negative, 1527.1 HIV positive, and 17.6 recently infected weighted cases [@MPHIA2016_FirstReport].
#' We use these national counts to generate simulated district-level counts of recent infections for the
#' purpose of demonstrating the small-area incidence model.
#'
#' The expected number of 'recent'
#'
##+ simulate data, evaluate=FALSE
#'

##+ show data, echo=FALSE
tab_data <- mw@data[c("district", "pop15to49", "adultart", "nsamp", "npos", "nrecent", "prev_est", "prev_se")]
tab_data$pop15to49 <- round(tab_data$pop15to49, -3)
tab_data$adultart <- round(tab_data$adultart, -2)
tab_data$prev_est <- round(100*tab_data$prev_est, 1)
tab_data$prev_se <- round(100*tab_data$prev_se, 1)
tab_data[["ANC prev"]] <- with(mw@data, round(100 * ancrt_pos / ancrt_n, 1))
tab_data[["ANC artcov"]] <- with(mw@data, round(100 * ancrt_art / ancrt_pos, 1))

knitr::kable(tab_data[distord,], row.names=FALSE, align="lccccccccc",
             caption="Simulated district-level data iputs")
             
