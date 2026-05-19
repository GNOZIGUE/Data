

# ── 0. Installation ───────────────────────────────────────────────────────────
pkgs <- c(
  "readxl", "dplyr", "tidyverse", "brms", "loo",
  "tidyr", "ggplot2", "sf", "bayesplot", "posterior",
  "viridis", "ggspatial",
  "cmdstanr",       # backend rapide (remplace rstan, gain x2–x5)
  "future",         # parallélisme
  "future.apply",
  "parallel"
)
new_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs)

# Installation CmdStan (à faire UNE seule fois)
library(cmdstanr)
if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  # Sur Windows : CmdStan nécessite RTools (https://cran.r-project.org/bin/windows/Rtools/)
  # Vérifier que RTools est installé avant cette ligne
  cmdstanr::install_cmdstan(cores = parallel::detectCores())
}

lapply(pkgs, library, character.only = TRUE)


# ── 1. Détection des ressources i7-8750H ──────────────────────────────────────
n_physical <- parallel::detectCores(logical = FALSE)  # = 6 cœurs physiques i7-8750H (fixe pour cette machine)
n_logical  <- parallel::detectCores(logical = TRUE)   # = 12 (HT activé)
n_chains   <- 4L

# Sur i7-8750H : 12 threads logiques / 4 chaînes = 3 threads/chaîne
# Stan utilisera les 12 threads HT efficacement avec reduce_sum
threads_per_chain <- floor(n_logical / n_chains)  # = 3

message(sprintf(
  " Cœurs physiques: %d | Logiques (HT): %d | Threads/chaîne: %d",
  n_physical, n_logical, threads_per_chain
))
# Affichera : Cœurs physiques: 6 | Logiques (HT): 12 | Threads/chaîne: 3


# ── 2. Plan future WINDOWS ────────────────────────────────────────────────────
# IMPORTANT : sur Windows, utiliser OBLIGATOIREMENT multisession (pas multicore)
# multicore utilise fork() qui n'existe pas sous Windows → erreur

# Pour les diagnostics (LOO/kfold) en parallèle APRÈS les fits
future::plan(future::multisession, workers = 2L)
# workers = 2 : un worker par modèle pour les diagnostics parallèles


# ── 3. Données ────────────────────────────────────────────────────────────────
precipitation <- readxl::read_excel(
  "../Matrice_des_parametres.xlsx",
  sheet = "Pr\u00E9cipitation"
)

precip <- precipitation %>%
  pivot_longer(
    cols      = starts_with("20"),
    names_to  = "Year",
    values_to = "precip_max"
  ) %>%
  mutate(
    year    = rep(seq(24, 1, -1), 1568),
    station = rep(1:1568, each = 24)
  ) %>%
  mutate(
    year_scaled = as.vector(scale(year)),
    LAT_scaled  = as.vector(scale(LAT)),
    LON_scaled  = as.vector(scale(LON))
  )


# ── 4. Formules et Priors ─────────────────────────────────────────────────────
formula_mixte <- bf(
  precip_max ~ 1 + LAT_scaled + LON_scaled + year_scaled + (1 | station),
  sigma      ~ 1 + LAT_scaled + LON_scaled + year_scaled + (1 | station),
  xi         ~ 1
)

# Option A : GP exact (lent)
formula_gp_exact <- bf(
  precip_max ~ 1 + LAT_scaled + LON_scaled + year_scaled +
    gp(LAT_scaled, LON_scaled, scale = FALSE) + (1 | station),
  sigma      ~ 1 + LAT_scaled + LON_scaled + year_scaled +
    gp(LAT_scaled, LON_scaled, scale = FALSE) + (1 | station),
  xi         ~ 1
)

# Option B : GP approché (RECOMMANDÉ)
# c = 5/4 active l'approximation de Hilbert (HSGP) :
# - Réduit la complexité de O(n³) à O(n log n)
# - Résultats quasi-identiques pour données spatiales régulières
#
formula_gp_approx <- bf(
  precip_max ~ 1 + LAT_scaled + LON_scaled + year_scaled +
    gp(LAT_scaled, LON_scaled, scale = FALSE, c = 5/4, k = 10) + (1 | station),
  sigma      ~ 1 + LAT_scaled + LON_scaled + year_scaled +
    gp(LAT_scaled, LON_scaled, scale = FALSE, c = 5/4, k = 10) + (1 | station),
  xi         ~ 1
)

# Choisir la formule GP à utiliser (changer ici si besoin)
formula_gp <- formula_gp_exact   

priors_communs <- c(
  prior(normal(49.44304, 53.74601), class = "Intercept"),
  prior(logistic(log(48.27761), 1), class = "Intercept", dpar = "sigma"),
  prior(normal(0, 0.5),            class = "Intercept", dpar = "xi"),
  prior(normal(0, 5),              class = "b"),
  prior(normal(0, 5),              class = "b",          dpar = "sigma"),
  prior(exponential(1),            class = "sd"),
  prior(exponential(1),            class = "sd",         dpar = "sigma")
)


# ── 5. Paramètres MCMC communs ────────────────────────────────────────────────
mcmc_args <- list(
  chains  = n_chains,          # 4 chaînes
  iter    = 2000L,
  warmup  = 1000L,
  cores   = n_chains,          # 4 cœurs pour les chaînes (1 par chaîne)
  threads = threading(threads_per_chain),   # 3 threads Stan par chaîne
  backend = "cmdstanr",        # OBLIGATOIRE pour threading + rapidité
  seed    = 1234L,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)
# Bilan d'utilisation CPU :
# 4 chaînes × 3 threads = 12 threads logiques → 100 % utilisation HT i7-8750H


# ── 6. Fit SÉQUENTIEL (recommandé sur 6 cœurs physiques) ─────────────────────
# Pourquoi séquentiel et non parallèle entre modèles ?
# → Lancer 2 modèles en même temps signifie 2 × 4 chaînes = 8 chaînes
#   pour seulement 6 cœurs physiques → contention, temps x2 sans gain.
# → En séquentiel, chaque modèle dispose de TOUS les 12 threads HT.

message("\n[Fit 1/2] Modèle Mixte GEV — démarrage...")
t1 <- proc.time()

fit_gev_mixte <- do.call(brm, c(
  list(formula = formula_mixte,
       data    = precip,
       family  = gen_extreme_value(),
       prior   = priors_communs),
  mcmc_args
))

t_mixte <- proc.time() - t1
message(sprintf("[Fit 1/2] Terminé en %.1f min", t_mixte[3] / 60))

# Sauvegarde intermédiaire (protection contre plantage/coupure courant)
save.image(file = "fit_gev_mixte.RData")


message("\n[Fit 2/2] Modèle GP GEV — démarrage...")
t2 <- proc.time()

fit_complete <- do.call(brm, c(
  list(formula = formula_gp,
       data    = precip,
       family  = gen_extreme_value(),
       prior   = priors_communs),
  mcmc_args
))

t_gp <- proc.time() - t2
message(sprintf("[Fit 2/2] Terminé en %.1f min", t_gp[3] / 60))

save.image(file = "fit_complete.RData")


# ── 7. Diagnostics PARALLÈLES (multisession Windows) ─────────────────────────
# Maintenant que les fits sont terminés, on peut lancer les diagnostics
# en parallèle : chaque worker traite UN modèle.
message("\n[Diagnostic] LOO / WAIC / K-fold en parallèle...")
t_diag <- proc.time()

loo_job <- future::future({
  list(
    mixte = loo(fit_gev_mixte, cores = 2L),
    gp    = loo(fit_complete,  cores = 2L)
  )
}, seed = TRUE)

waic_job <- future::future({
  list(
    mixte = waic(fit_gev_mixte),
    gp    = waic(fit_complete)
  )
}, seed = TRUE)

kfold_job <- future::future({
  list(
    mixte = kfold(fit_gev_mixte, K = 10, save_fits = FALSE, cores = 2L),
    gp    = kfold(fit_complete,  K = 10, save_fits = FALSE, cores = 2L)
  )
}, seed = TRUE)

loo_res   <- future::value(loo_job)
waic_res  <- future::value(waic_job)
kfold_res <- future::value(kfold_job)

loo_mixte  <- loo_res$mixte;   loo_gp   <- loo_res$gp
waic_mixte <- waic_res$mixte;  waic_gp <- waic_res$gp
kfold_mixte   <- kfold_res$mixte; kfold_gp    <- kfold_res$gp

t_diag_end <- proc.time() - t_diag
message(sprintf("[Diagnostic] Terminé en %.1f min", t_diag_end[3] / 60))


# ── 8. Comparaison & Résumé ───────────────────────────────────────────────────
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


# ── 9. Nettoyage et sauvegarde ────────────────────────────────────────────────
future::plan(future::sequential)   # libère les workers multisession
save.image(file = "Kossivi_BHM_windows.RData")
message("[Done] Tout sauvegardé dans Kossivi_BHM_windows.RData")
