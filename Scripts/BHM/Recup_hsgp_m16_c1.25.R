
library(readxl); library(brms);library(cmdstanr); library(ggplot2)
library(brms); library(cmdstanr);library(tidyverse); library(sf)
library(bayesplot);library(posterior);library(loo);library(viridis) ; library(ggspatial);library(tidyr);library(dplyr)
library(future)

load("fit_hsgp_m16_c125.rds")

fit1 <- readRDS("C:/Users/DELL/Desktop/Parallelisme/Parallel_HSGP/Result_Kossivi_BHM_HSGP_GP_Serveur/fit_mod_chain_1.rds")
fit2 <- readRDS("C:/Users/DELL/Desktop/Parallelisme/Parallel_HSGP/Result_Kossivi_BHM_HSGP_GP_Serveur/fit_mod_chain_2.rds")
fit3 <- readRDS("C:/Users/DELL/Desktop/Parallelisme/Parallel_HSGP/Result_Kossivi_BHM_HSGP_GP_Serveur/fit_mod_chain_3.rds")
fit4 <- readRDS("C:/Users/DELL/Desktop/Parallelisme/Parallel_HSGP/Result_Kossivi_BHM_HSGP_GP_Serveur/fit_mod_chain_4.rds")

fit_hsgp <- combine_models(fit1, fit2, fit3, fit4, check_data = TRUE)

# Doit afficher : chains = 4
print(fit_hsgp)

# Vérification diagnostique complète
summary(fit_hsgp)

# R-hat et ESS (tout R-hat doit être < 1.01)
posterior::summarise_draws(
  as_draws(fit_hsgp),
  "mean", "sd", "rhat"
)

# Trace plots pour visualiser la convergence inter-chaînes
plot(fit_hsgp, ask = FALSE)

shinystan::launch_shinystan(fit_hsgp)  # traceplots 
plot(conditional_effects(fit_hsgp))  # non-stationnarité spatiale

# Si ShinyStan ne le détecte pas automatiquement, préparez :
y_obs <- fit_hsgp$data$precip_max
# Extraction standard — toutes les observations
yrep_m16 <- posterior_predict(fit_hsgp)

""" 
% latex table generated in R 4.5.1 by xtable 1.8-4 package
% Tue May 19 00:01:05 2026
\begin{table}[ht]
\centering
\begin{tabular}{lrrrrrrr}
\toprule
Parameter & Rhat & n\_eff & mean & sd & se\_mean & 50\% & 97.5\% \\ 
\midrule
b\_Intercept & 1.001 & 2852 & 28.591 & 4.000 & 0.075 & 28.687 & 36.371 \\ 
b\_sigma\_Intercept & 1.001 & 2052 & 2.649 & 0.160 & 0.004 & 2.644 & 2.978 \\ 
b\_xi\_Intercept & 1.000 & 4000 & 0.167 & 0.004 & 0.000 & 0.167 & 0.174 \\ 
b\_LAT\_scaled & 1.001 & 3133 & -2.740 & 2.685 & 0.048 & -2.790 & 2.574 \\ 
b\_LON\_scaled & 1.001 & 2990 & -0.695 & 2.628 & 0.048 & -0.672 & 4.512 \\ 
b\_year\_scaled & 0.999 & 4000 & 6.118 & 0.077 & 0.001 & 6.117 & 6.268 \\ 
b\_sigma\_LAT\_scaled & 1.001 & 1343 & 0.050 & 0.146 & 0.004 & 0.035 & 0.385 \\ 
b\_sigma\_LON\_scaled & 1.000 & 2500 & -0.181 & 0.130 & 0.003 & -0.178 & 0.078 \\ 
b\_sigma\_year\_scaled & 0.999 & 4000 & 0.427 & 0.005 & 0.000 & 0.427 & 0.438 \\ 
sd\_station\_\_Intercept & 1.007 & 1254 & 0.085 & 0.065 & 0.002 & 0.072 & 0.243 \\ 
sd\_station\_\_sigma\_Intercept & 1.001 & 1950 & 0.006 & 0.004 & 0.000 & 0.005 & 0.017 \\ 
sdgp\_gpLAT\_scaledLON\_scaled & 1.001 & 2054 & 56.633 & 15.601 & 0.344 & 53.225 & 95.229 \\ 
sdgp\_sigma\_gpLAT\_scaledLON\_scaled & 1.001 & 1796 & 1.709 & 0.673 & 0.016 & 1.515 & 3.580 \\ 
lscale\_gpLAT\_scaledLON\_scaled & 1.001 & 3044 & 0.132 & 0.039 & 0.001 & 0.129 & 0.215 \\ 
lscale\_sigma\_gpLAT\_scaledLON\_scaled & 1.002 & 1310 & 0.180 & 0.081 & 0.002 & 0.167 & 0.364 \\ 
Intercept & 1.001 & 2852 & 28.591 & 4.000 & 0.075 & 28.687 & 36.371 \\ 
Intercept\_sigma & 1.001 & 2052 & 2.649 & 0.160 & 0.004 & 2.644 & 2.978 \\ 
Intercept\_xi & 1.000 & 4000 & 0.167 & 0.004 & 0.000 & 0.167 & 0.174 \\ 
\bottomrule
\end{tabular}
\end{table}
"""

# ── 1. Extraire les draws sous forme de tableau 3D (iterations × chains × variables) ──
draws_hsgp <- as_draws_array(fit_hsgp)

# ── 2. Table de correspondance : nom brms → expression R (syntaxe plotmath) ──────────
param_labels <- c(
  # Paramètre μ (location)
  "b_Intercept"                    = "beta[0]^mu",
  "b_LAT_scaled"                   = "beta[1]^mu",
  "b_LON_scaled"                   = "beta[2]^mu",
  "b_year_scaled"                  = "beta[3]^mu",
  
  # Paramètre σ (scale, log-link)
  "b_sigma_Intercept"              = "beta[0]^sigma",
  "b_sigma_LAT_scaled"             = "beta[1]^sigma",
  "b_sigma_LON_scaled"             = "beta[2]^sigma",
  "b_sigma_year_scaled"            = "beta[3]^sigma",
  
  # Paramètre ξ (shape)
  "b_xi_Intercept"                 = "beta[0]^xi",
  
  # Effets aléatoires (écarts-types inter-station)
  "sd_station__Intercept"          = "u^mu",
  "sd_station__sigma_Intercept"    = "u^sigma",
  
  # Hyperparamètres du GP (si vous voulez les inclure)
  "sdgp_gp_LAT_scaled_LON_scaled"              = "alpha[GP]^mu",
  "lscale_gp_LAT_scaled_LON_scaled"            = "ell^mu",
  "sdgp_sigma_gp_LAT_scaled_LON_scaled"        = "alpha[GP]^sigma",
  "lscale_sigma_gp_LAT_scaled_LON_scaled"      = "ell^sigma"
)

# ── 3. Renommer les dimensions du tableau draws ───────────────────────────────────────
old_names <- dimnames(draws_hsgp)$variable
dimnames(draws_hsgp)$variable <- ifelse(
  old_names %in% names(param_labels),
  param_labels[old_names],
  old_names   # conserver le nom original si absent du mapping
)


# ── 4. Sélectionner uniquement les paramètres fixes (facultatif, pour lisibilité) ─────
fixed_pars <- c(
  "beta[0]^mu", "beta[1]^mu", "beta[2]^mu", "beta[3]^mu",
  "beta[0]^sigma", "beta[1]^sigma", "beta[2]^sigma", "beta[3]^sigma",
  "beta[0]^xi", "u^mu", "u^sigma"
)

draws_fixed_hsgp <- subset_draws(draws_hsgp, variable = fixed_pars)



fix_labels <- function(p) {
  p$facet$params$labeller <- label_parsed
  p
}
# Densités postérieures
fix_labels(
  mcmc_dens_overlay(draws_fixed_hsgp) + 
    ggtitle("Posterior densities of parameters") +
    theme(plot.title = element_text(hjust = 0.5))
)

# Traces MCMC
fix_labels(
  mcmc_trace(draws_fixed_hsgp) + 
    ggtitle("MCMC traces of parameters") +
    theme(plot.title = element_text(hjust = 0.5))
)

# ACF
fix_labels(
  mcmc_acf(draws_fixed_hsgp) + 
    ggtitle("Autocorrélation (ACF)") +
    theme(plot.title = element_text(hjust = 0.5))
)



#color_scheme_set("blue")  # Tous bleus
#Rhat
mcmc_plot(fit_hsgp, type = "rhat_hist")+
  ggtitle("rhat histogram") +
  theme(plot.title = element_text(hjust = 0.5))

mcmc_plot(fit_hsgp, type = "rhat")
#acf
mcmc_plot(fit_hsgp, type = "acf")+
  ggtitle("acf") +
  theme(plot.title = element_text(hjust = 0.1))
#mcmc_plot(fit_hsgp, type = "acf_bar")

# Differentes palettes de couleurs
color_scheme_set("mix-blue-red") 
color_scheme_set("brewer-Spectral")
# Définir manuellement 6 couleurs (bayesplot en utilise 6 niveaux)
color_scheme_set(c(
  "#E69F00",   # light  (chaîne 1 claire)
  "#D55E00",   # light mid
  "#CC79A7",   # mid
  "#0072B2",   # mid dark
  "#009E73",   # dark mid
  "#56B4E9"    # dark   (chaîne 4 foncée)
))

mcmc_plot(fit_hsgp, type = "dens_overlay") +
  ggtitle("Posterior densities of coefficients") +
  theme(plot.title = element_text(hjust = 0.5))

mcmc_plot(fit_hsgp, type = "dens_overlay") +
  scale_color_manual(
    values = c(
      "chain:1" = "#E41A1C",
      "chain:2" = "#377EB8", 
      "chain:3" = "#4DAF4A",
      "chain:4" = "#984EA3"
    )
  ) +
  ggtitle("Posterior densities of coefficients") +
  theme(plot.title = element_text(hjust = 0.5))



mcmc_plot(fit_hsgp, type = "dens_overlay") +
  ggtitle("Posterior densities of coefficients") +
  theme(plot.title = element_text(hjust = 0.5))

mcmc_plot(fit_hsgp, type = "trace")+
  ggtitle("Posterior traces of parameters") +
  theme(plot.title = element_text(hjust = 0.5))



## 4.4 Dist Résidus---------
residuals_hsgp <- residuals(fit_hsgp, summary = FALSE)
residuals_mean <- apply(residuals_hsgp, 2, mean)

ggplot(data.frame(residual = residuals_mean)) +
  geom_histogram(aes(x = residual), bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Distribution of mean Bayesian residuals", 
       x = "Average residue", y = "Frequency") +
  theme_minimal()


# 5. validation prédictive postérieure --------
# 5.1 PPC pour les statistiques principales
pp_check(fit_hsgp, type = "dens_overlay", ndraws = 100) +
  ggtitle("Posterior Predictive Check: Distribution")  +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

# Identifier les valeurs anormalement élevées
precip_20max = precip %>%
  arrange(desc(precip_max)) %>%
  head(20) %>%
  select(STATIONS, Year, precip_max)
# Des valeurs > 500 mm méritent vérification de l'unité ou de la saisie
# Distribution par station
precip %>%
  group_by(STATIONS) %>%
  summarise(max_obs = max(precip_max), n = n()) %>%
  arrange(desc(max_obs))
# Quelques stations "hyperpluvieuses" tirent-elles toute la queue ?



pp_check(fit_hsgp, type = "stat", stat = "mean", ndraws = 100) +
  ggtitle("Posterior Predictive Check: Average")

pp_check(fit_hsgp, type = "stat", stat = "sd", ndraws = 100) +
  ggtitle("Posterior Predictive Check: Standard deviation")

# 5.2 PPC pour les extrêmes (importante pour GEV)
pp_check(fit_hsgp, type = "stat", stat = function(y) quantile(y, 0.95)) +
  ggtitle("Posterior Predictive Check: 95e quantile")

pp_check(fit_hsgp, type = "stat", stat = function(y) max(y)) +
  ggtitle("Posterior Predictive Check: Maximum")

loo_m16 <-  loo(fit_hsgp)
waic_m16 <- waic(fit_hsgp)



# 7.Predict-----

# newdata
newdata <- precip %>% 
  distinct(STATIONS, LAT, LON, LAT_scaled, LON_scaled, year_scaled, Year, station)


# Posibilité 1: Tirages postérieurs de mu, sigma, xi pour chaque horizon
get_gev_params <- function(fit_hsgp, newdata) {
  mu    <- posterior_epred(fit_hsgp, newdata = newdata, dpar = "mu",    re_formula = NULL)
  sigma <- posterior_epred(fit_hsgp, newdata = newdata, dpar = "sigma", re_formula = NULL)
  xi    <- posterior_epred(fit_hsgp, newdata = newdata, dpar = "xi",    re_formula = NULL)
  list(mu = mu, sigma = sigma, xi = xi)
}

params  <- get_gev_params(fit_hsgp, newdata)


mu_med   <- apply(params$mu, 2, median)
mu_q25   <- apply(params$mu, 2, quantile, 0.025)
mu_q97_5  <- apply(params$mu, 2, quantile, 0.975)

sigma_med   <- apply(params$sigma, 2, median)
sigma_q25   <- apply(params$sigma, 2, quantile, 0.025)
sigma_q97_5  <- apply(params$sigma, 2, quantile, 0.975)

xi_med   <- apply(params$xi, 2, median)
xi_q25   <- apply(params$xi, 2, quantile, 0.025)
xi_q97_5  <- apply(params$xi, 2, quantile, 0.975)

# 8. RL ----------
# Cas 1
gev_return_level <- function(mu, sigma, xi, T) {
  # Probabilité de non-dépassement annuelle
  p  <- 1 - 1/T
  yp <- -log(p)           
  
  # Formule GEV (xi != 0)
  ifelse(
    abs(xi) < 1e-6,
    mu - sigma * log(yp),                          # cas limite Gumbel
    mu + sigma / xi * (yp^(-xi) - 1)              # cas général GEV
  )
}

# Calcul sur chaque tirage
rl_50  <- gev_return_level(mu_med, sigma_med, xi_med, T = 50)
rl_50_q025  <- gev_return_level(mu_q25, sigma_q25, xi_q25, T = 50)
rl_50_q975  <- gev_return_level(mu_q97_5, sigma_q97_5, xi_q97_5, T = 50)

rl_100 <- gev_return_level(mu_med, sigma_med, xi_med, T = 100)

summarise_rl <- function(rl, T) {
  tibble(
    period   = paste0(T, " years"),
    min      = min(rl),
    q025     = quantile(rl, 0.025),
    median   = median(rl),
    mean     = mean(rl),
    q975     = quantile(rl, 0.975),
    max      = max(rl)
  )
}

bind_rows(
  summarise_rl(rl_50,  50),
  summarise_rl(rl_100, 100)
)

# 9. Carte cas1 -------- 
res_rl = data.frame(
  LAT = precip$LAT,
  LON = precip$LON,
  Year = precip$Year,
  rl_50 = rl_50,
  rl_100 = rl_100,
  rl_50_q025 = rl_50_q025,
  rl_50_q975 = rl_50_q975
)

ggplot() +
  # Couleur des frontières
  geom_point(data = res_rl, aes(x = LON, y = LAT, color =  rl_50),
             size = 5, stat = "identity", shape = "square") + 
  scale_color_gradient(low = c("lightblue", "darkblue", "green"), high = c("yellow", "orange", "red")) + # Dégradé de couleur
  
  #scale_color_gradientn(colors = c("lightblue", "darkblue", "green", "yellow", "orange", "red")) + # Dégradé de couleur
  #scale_fill_viridis_d(option = "A")+
  
  scale_x_continuous(name = "Longitude", 
                     breaks = seq(-17, 11, by = 12),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(name = "Latitude", 
                     breaks = seq(4, 20, by = 10),
                     labels = function(y) paste0(y, "°N")) +
  
  labs(title = "50-year return level", 
       #  x = "Longitude",
       #    y = "Latitude",
       color = "mm/jr")+
  theme_minimal() +
  facet_wrap(~ Year) +
  
  theme(legend.title = element_text(size = 8, face = "bold", ), 
        legend.position="right", legend.text = element_text(size = 10, face = "bold"),
        #legend.position.inside =c(15,18), 
        plot.title = element_text(hjust = 0.5, size = 10, ), 
        strip.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"), 
        axis.title.y = element_text(size = 10, face = "bold"))+ 
  
  coord_fixed(xlim = c(-16.75, 10.1), ylim = c(4.25, 17.75)) 

ggplot() +
  # Couleur des frontières
  geom_point(data = res_rl, aes(x = LON, y = LAT, color =  rl_100),
             size = 5, stat = "identity", shape = "square") + 
  scale_color_gradient(low = c("lightblue", "darkblue", "green"), high = c("yellow", "orange", "red")) + # Dégradé de couleur
  
  #scale_color_gradientn(colors = c("lightblue", "darkblue", "green", "yellow", "orange", "red")) + # Dégradé de couleur
  #scale_fill_viridis_d(option = "A")+
  
  scale_x_continuous(name = "Longitude", 
                     breaks = seq(-17, 11, by = 12),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(name = "Latitude", 
                     breaks = seq(4, 20, by = 10),
                     labels = function(y) paste0(y, "°N")) +
  
  labs(title = "100-years return level", 
       #  x = "Longitude",
       #    y = "Latitude",
       color = "mm/jr")+
  theme_minimal() +
  facet_wrap(~ Year) +
  
  theme(legend.title = element_text(size = 8, face = "bold", ), 
        legend.position="right", legend.text = element_text(size = 10, face = "bold"),
        #legend.position.inside =c(15,18), 
        plot.title = element_text(hjust = 0.5, size = 10, ), 
        strip.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"), 
        axis.title.y = element_text(size = 10, face = "bold"))+ 
  
  coord_fixed(xlim = c(-16.75, 10.1), ylim = c(4.25, 17.75)) 


borders("world", regions = c("Benin", "Burkina Faso", "Cabo Verde", "Ivory Coast", 
                             "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
                             "Liberia", "Mali", "Mauritania", "Niger", "Nigeria", 
                             "Senegal", "Sierra Leone", "Togo"),
        colour = "black",size=.8) 



#

save.image("fit_hsgp_m16_c125.rds")



