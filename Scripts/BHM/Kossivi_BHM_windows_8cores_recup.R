#1. Packages---------
library(readxl); library(brms);library(cmdstanr); library(ggplot2)
library(brms); library(cmdstanr);library(tidyverse); library(sf)
library(bayesplot);library(posterior);library(loo);library(viridis) ; library(ggspatial);library(tidyr);library(dplyr)
library(future)


#2. Chargement de l'espace de travail------
load("Kossivi_BHM_8cores_final.RData") 

#Récupération des données des modèles
load("fit_complete.RData")
load("fit_gev_mixte.RData")
fit_gev_mixte
fit_complete

#3. Diagnostic bayésien -------------
t_diag <- proc.time()

loo_job <- list(
    mixte = loo(fit_gev_mixte),
    gp    = loo(fit_complete)
  )


waic_job <- list(
    mixte = waic(fit_gev_mixte),
    gp    = waic(fit_complete)
  )

# Attention très lent
kfold_job <- list(
    mixte = kfold(fit_gev_mixte, K = 10),
    gp_kfold    = kfold(fit_complete,  K = 10)
  )

#Recupération des resultats
loo_res   <- future::value(loo_job)
waic_res  <- future::value(waic_job)
kfold_res <- future::value(kfold_job)

loo_mixte  <- loo_res$mixte;   loo_gp   <- loo_res$gp
waic_mixte <- waic_res$mixte;  waic_gp <- waic_res$gp
kfold_mixte   <- kfold_res$mixte; kfold_gp    <- kfold_res$gp

t_diag_end <- proc.time() - t_diag
message(sprintf("[Diagnostic] Terminé en %.1f min", t_diag_end[3] / 60))


# Comparaison & Résumé 
cat("\n=== Comparaison LOO ===\n")
print(loo_compare(loo_mixte, loo_gp))

t_total <- t_mixte[3] + t_gp[3] + t_diag_end[3]
message(sprintf(
  "\n[Bilan] Temps total : %.1f h (Mixte: %.1f min | GP: %.1f min | Diag: %.1f min)",
  t_total / 3600,
  t_mixte[3] / 60,
  t_gp[3] / 60,
  t_diag_end[3] / 60
))




#4. Trace plots et diagnostics MCMC---------
#mcmc_plots <- plot(fit_complete, ask = FALSE)


## 4.1 R-hat et ESS (Effective Sample Size)-------
rhats_m <- rhat(fit_gev_mixte)
ess_bulk <- ess_bulk(fit_gev_mixte)
ess_tail <- ess_tail(fit_gev_mixte)


rhats_GP <- rhat(fit_complete)
ess_bulk_GP <- ess_bulk(fit_complete)
ess_tail_GP <- ess_tail(fit_complete)

cat("=== Diagnostics MCMC ===\n")
cat(sprintf("R-hat max: %.3f\n", max(rhats_m, na.rm = TRUE)))
cat(sprintf("R-hat > 1.05: %d\n", sum(rhats_m > 1.05, na.rm = TRUE)))
cat(sprintf("ESS bulk min: %.0f\n", min(ess_bulk, na.rm = TRUE)))
cat(sprintf("ESS tail min: %.0f\n", min(ess_tail, na.rm = TRUE)))

cat(sprintf("R-hat max: %.3f\n", max(rhats_GP, na.rm = TRUE)))
cat(sprintf("R-hat > 1.05: %d\n", sum(rhats_GP > 1.05, na.rm = TRUE)))
cat(sprintf("ESS bulk min: %.0f\n", min(ess_bulk_GP, na.rm = TRUE)))
cat(sprintf("ESS tail min: %.0f\n", min(ess_tail_GP, na.rm = TRUE)))

##4.2 plots-----



# ── 1. Extraire les draws sous forme de tableau 3D (iterations × chains × variables) ──
draws <- as_draws_array(fit_complete)

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
old_names <- dimnames(draws)$variable
dimnames(draws)$variable <- ifelse(
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

draws_fixed <- subset_draws(draws, variable = fixed_pars)



fix_labels <- function(p) {
  p$facet$params$labeller <- label_parsed
  p
}
# Densités postérieures
fix_labels(
  mcmc_dens_overlay(draws_fixed) + 
    ggtitle("Posterior densities of parameters") +
    theme(plot.title = element_text(hjust = 0.5))
)

# Traces MCMC
fix_labels(
  mcmc_trace(draws_fixed) + 
    ggtitle("MCMC traces of parameters") +
    theme(plot.title = element_text(hjust = 0.5))
)

# ACF
fix_labels(
  mcmc_acf(draws_fixed) + 
    ggtitle("Autocorrélation (ACF)") +
    theme(plot.title = element_text(hjust = 0.5))
)



#color_scheme_set("blue")  # Tous bleus
#Rhat
mcmc_plot(fit_complete, type = "rhat_hist")+
  ggtitle("rhat histogram") +
  theme(plot.title = element_text(hjust = 0.5))

mcmc_plot(fit_complete, type = "rhat")
#acf
mcmc_plot(fit_complete, type = "acf")+
  ggtitle("acf") +
  theme(plot.title = element_text(hjust = 0.1))
#mcmc_plot(fit_complete, type = "acf_bar")

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

mcmc_plot(fit_complete, type = "dens_overlay") +
  ggtitle("Posterior densities of coefficients") +
  theme(plot.title = element_text(hjust = 0.5))

mcmc_plot(fit_complete, type = "dens_overlay") +
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



mcmc_plot(fit_complete, type = "dens_overlay") +
  ggtitle("Posterior densities of coefficients") +
  theme(plot.title = element_text(hjust = 0.5))

mcmc_plot(fit_complete, type = "trace")+
  ggtitle("Posterior traces of parameters") +
  theme(plot.title = element_text(hjust = 0.5))



## 4.4 Dist Résidus---------
residuals <- residuals(fit_complete, summary = FALSE)
residuals_mean <- apply(residuals, 2, mean)

ggplot(data.frame(residual = residuals_mean)) +
  geom_histogram(aes(x = residual), bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Distribution of mean Bayesian residuals", 
       x = "Average residue", y = "Frequency") +
  theme_minimal()


# 5. validation prédictive postérieure --------
# 5.1 PPC pour les statistiques principales
pp_check(fit_complete, type = "dens_overlay", ndraws = 100) +
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



pp_check(fit_complete, type = "stat", stat = "mean", ndraws = 100) +
  ggtitle("Posterior Predictive Check: Average")

pp_check(fit_complete, type = "stat", stat = "sd", ndraws = 100) +
  ggtitle("Posterior Predictive Check: Standard deviation")

# 5.2 PPC pour les extrêmes (importante pour GEV)
pp_check(fit_complete, type = "stat", stat = function(y) quantile(y, 0.95)) +
  ggtitle("Posterior Predictive Check: 95e quantile")

pp_check(fit_complete, type = "stat", stat = function(y) max(y)) +
  ggtitle("Posterior Predictive Check: Maximum")


#6. SHINY -------------
shinystan::launch_shinystan(fit_complete)  # traceplots 
plot(conditional_effects(fit_complete))  # non-stationnarité spatiale

# 7.Predict-----

# newdata
newdata <- precip %>% 
  distinct(STATIONS, LAT, LON, LAT_scaled, LON_scaled, year_scaled, Year, station)


# Posibilité 1: Tirages postérieurs de mu, sigma, xi pour chaque horizon
get_gev_params <- function(fit, newdata) {
  mu    <- posterior_epred(fit, newdata = newdata, dpar = "mu",    re_formula = NULL)
  sigma <- posterior_epred(fit, newdata = newdata, dpar = "sigma", re_formula = NULL)
  xi    <- posterior_epred(fit, newdata = newdata, dpar = "xi",    re_formula = NULL)
  list(mu = mu, sigma = sigma, xi = xi)
}

params  <- get_gev_params(fit_complete, newdata)


mu_med   <- apply(params$mu, 2, median)
mu_q25   <- apply(params$mu, 2, quantile, 0.025)
mu_q97_5  <- apply(params$mu, 2, quantile, 0.975)

sigma_med   <- apply(params$sigma, 2, median)
sigma_q25   <- apply(params$sigma, 2, quantile, 0.025)
sigma_q97_5  <- apply(params$sigma, 2, quantile, 0.975)

xi_med   <- apply(params$xi, 2, median)
xi_q25   <- apply(params$xi, 2, quantile, 0.025)
xi_q97_5  <- apply(params$xi, 2, quantile, 0.975)


#Possibilité2 : Tirages MCMC pour chaque paramètre séparément
mu_draws <- fitted(
  fit_complete,
  newdata = newdata,
  dpar = "mu",
  re_formula = NULL,
  allow_new_levels = TRUE,
  summary = FALSE
)

sigma_draws <- fitted(
  fit_complete,
  newdata = newdata,
  dpar = "sigma",
  re_formula = NULL,
  allow_new_levels = TRUE,
  summary = FALSE
)

xi_draws <- fitted(
  fit_complete,
  newdata = newdata,
  dpar = "xi",
  re_formula = NULL,
  allow_new_levels = TRUE,
  summary = FALSE
)

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

#Cas 2
gev_return_level2 <- function(mu, sigma, xi, p) {
  # vectorisé
  eps <- 1e-6
  out <- ifelse(
    abs(xi) > eps,
    mu + (sigma/xi) * ((-log(p))^(-xi) - 1),
    mu - sigma * log(-log(p))
  )
  out
}

p50 <- 1 - 1/50
p100 <- 1 - 1/100

# pour chaque station (colonne) et chaque tirage (ligne)
zT_draws50 <- gev_return_level2(mu_draws, sigma_draws, xi_draws, p50)
zT_draws100 <- gev_return_level2(mu_draws, sigma_draws, xi_draws, p100)

# résumer par station: médiane et intervalle de crédibilité
zT_med50   <- apply(zT_draws50, 2, median)
zT_low50   <- apply(zT_draws50, 2, quantile, 0.025)
zT_high50  <- apply(zT_draws50, 2, quantile, 0.975)

zT_med100   <- apply(zT_draws100, 2, median)
zT_low100   <- apply(zT_draws100, 2, quantile, 0.025)
zT_high100  <- apply(zT_draws100, 2, quantile, 0.975)

res_return <- cbind(
  newdata,
  zT_med50 = zT_med50,
  zT_low50 = zT_low50,
  zT_high50 = zT_high50,
  zT_med100 = zT_med100,
  zT_low100 = zT_low100,
  zT_high100 = zT_high100
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
  
  labs(title = "50-years return level", 
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



# 10. Carte cas2 ---------
## 10.1. Valeur médiane --------
### 10.1.1. Sans IC --------
ggplot() +
  # Couleur des frontières
  geom_point(data = res_return, aes(x = LON, y = LAT, color =  zT_med50),
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
  
  labs(title = "50-years return level", 
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





ggplot() +
  # Couleur des frontières
  geom_point(data = res_return, aes(x = LON, y = LAT, color =  zT_med100),
             size = 5, stat = "identity", shape = "square") + 
  
  scale_color_gradient(low = c("lightblue", "darkblue", "green"), high = c("yellow", "orange", "red")) + # Dégradé de couleur
  
  #scale_color_gradientn(colors = c("lightblue","darkblue", "green", "yellow", "orange", "red")) + # Dégradé de couleur
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
  
  theme(legend.title = element_text(size = 10, face = "bold", ), 
        legend.position="right", legend.text = element_text(size = 10, face = "bold"),
        #legend.position.inside =c(15,18), 
        plot.title = element_text(hjust = 0.5, size = 12, ), 
        strip.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"))+ 
  
  coord_fixed( xlim = c(-16.75, 10.1), ylim = c(4.25, 17.75))


###10.1.2 Avec IC ------


ggplot() +
  # 1. COULEUR = Médiane (votre carte existante)
  geom_point(data = res_return, 
             aes(x = LON, y = LAT, color = zT_med50),
             size = 5, shape = "square") +
  
  # 2. TAILLE = Incertitude (étendue IC 95%)
  scale_size_continuous(aes(size = zT_high50 - zT_low50), 
                        range = c(3, 12), name = "IC 95%") +
  
  # 3. BORDURES = Précision (inverse écart-type)
  geom_point(data = res_return, 
             aes(x = LON, y = LAT, color = zT_med50, size = 1/sd(zT_med50) ),
             shape = 1, stroke = 1.5, alpha = 0.8) +
  
  scale_color_gradient(low = c("lightblue", "darkblue", "green"), high = c("yellow", "orange", "red"), 
                       #midpoint = median(res_return$zT_med50),
                       name = "mm/day") +
  
  # Vos scales/échelles existants
  scale_x_continuous(name = "Longitude", breaks = seq(-17, 11, by = 12),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) +
  scale_y_continuous(name = "Latitude", breaks = seq(4, 20, by = 10),
                     labels = function(y) paste0(y, "°N")) +
  
  labs(title = "50-years return levels with Bayesian uncertainties",
       subtitle = "Colour = Median | Size = 95% CI width | Border = 1/sd") +
  
  theme_minimal() + facet_wrap(~ Year) +
  
  theme(legend.title = element_text(size = 8, face = "bold", ), 
        legend.position="right", legend.text = element_text(size = 10, face = "bold"),
        #legend.position.inside =c(15,18), 
        plot.title = element_text(hjust = 0.5, size = 10, ), 
        strip.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"), 
        axis.title.y = element_text(size = 10, face = "bold"))+ 
  
  coord_fixed( xlim = c(-16.75, 10.1), ylim = c(4.25, 17.75))


ggplot() +
  # 1. COULEUR = Médiane (votre carte existante)
  geom_point(data = res_return, 
             aes(x = LON, y = LAT, color = zT_med100),
             size = 5, shape = "square") +
  
  # 2. TAILLE = Incertitude (étendue IC 95%)
  scale_size_continuous(aes(size = zT_high100 - zT_low100), 
                        range = c(3, 12), name = "IC 95% (mm/jr)") +
  
  # 3. BORDURES = Précision (inverse écart-type)
  geom_point(data = res_return, 
             aes(x = LON, y = LAT, color = zT_med100, size = 1/sd(zT_med100) ),
             shape = 1, stroke = 1.5, alpha = 0.8) +
  
  scale_color_gradient(low = c("lightblue", "darkblue", "green"), high = c("yellow", "orange", "red"), 
                       #midpoint = median(res_return$zT_med100),
                       name = "mm/day") +
  
  # Vos scales/échelles existants
  scale_x_continuous(name = "Longitude", breaks = seq(-17, 11, by = 6),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) +
  scale_y_continuous(name = "Latitude", breaks = seq(4, 20, by = 3),
                     labels = function(y) paste0(y, "°N")) +
  
  labs(title = "100-years return levels with Bayesian uncertainties",
       subtitle = "Colour = Median | Size = 95% CI width | Border = 1/sd") +
  
  theme_minimal() + facet_wrap(~ Year) +
  
  theme(legend.title = element_text(size = 10, face = "bold", ), 
        legend.position="right", legend.text = element_text(size = 10, face = "bold"),
        #legend.position.inside =c(15,18), 
        plot.title = element_text(hjust = 0.5, size = 12, ), 
        strip.text = element_text(size = 8, face = "bold"), 
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"))+ 
  
  coord_fixed( xlim = c(-16.75, 10.1), ylim = c(4.25, 17.75))




save.image(file = "Kossivi_BHM_8cores_final.RData")
