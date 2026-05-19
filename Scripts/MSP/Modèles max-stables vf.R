
load("maxstable")

# Données et packages------------
library(SpatialExtremes); library(readxl); library(tidyr); library(dplyr); library(RColorBrewer); library(ggplot2); library(latex2exp); library(corrplot)

PL <- read.delim("../AO.txt")
#View(PL)
attach(PL)

coord <- matrix(c(LAT,LON), ncol=2) ;  colnames(coord) <- c("lat", "lon")

D_Mat = matrix(c(A2000, A2001, A2002, A2003, A2004, A2005, A2006, A2007, A2008, A2009,
                 A2010, A2011, A2012, A2013, A2014 ,A2015 ,A2016, A2017, A2018, A2019, A2020, A2021, A2022 , A2023),
               ncol=1568, byrow=TRUE )
#View(D_Mat)

#1. Stations supposées indépendantes ---------
# Modèle spatial en supposant que les stations météorologiques sont mutuellement indépendantes 
# Cela revient à utiliser une vraisemblance composite spéciale connue sous le nom de vraisemblance d’indépendanceVarin et al. (2011)


### a) Modèle 1 --------
loc.form1 <- ~ lat + lon
scale.form1 <- ~ lat + lon
shape.form1 <- ~ 1

m1 = fitspatgev(D_Mat, coord, loc.form1, scale.form1, shape.form1, corr =TRUE)
#m1 <- fitspatgev(t(PL[,-c(1,2)]), coord, loc.form0, scale.form0, shape.form0,corr =TRUE)

m1_loc_intercept = m1$param[1]
m1_loc_lat = m1$param[2]
m1_loc_lon = m1$param[3]

m1_scale_intercept = m1$param[4]
m1_scale_lat = m1$param[5]
m1_scale_lon = m1$param[6]

m1_shape_intercept = m1$param[7]


### b) Modèle 2 --------
loc.form2 <- ~ lat + lon
scale.form2 <- ~ lat 
shape.form2 <- ~ 1

m2 = fitspatgev(D_Mat, coord, loc.form2, scale.form2, shape.form2, corr =TRUE)
#m2 <- fitspatgev(t(PL[,-c(1,2)]), coord, loc.form1, scale.form1, shape.form1)

m2_loc_intercept = m2$param[1]
m2_loc_lat = m2$param[2]
m2_loc_lon = m2$param[3]

m2_scale_intercept = m2$param[4]
m2_scale_lat = m2$param[5]

m2_shape_intercept = m2$param[6]

### c) Modèle 3 --------
loc.form3 <- ~ lat + lon
scale.form3 <- ~ lat + lon
shape.form3 <- ~ lat + lon

m3 = fitspatgev(D_Mat, coord, loc.form3, scale.form3, shape.form3, corr =TRUE)
#m1 <- fitspatgev(t(PL[,-c(1,2)]), coord, loc.form0, scale.form0, shape.form0,corr =TRUE)

m3_loc_intercept = m3$param[1]
m3_loc_lat = m3$param[2]
m3_loc_lon = m3$param[3]

m3_scale_intercept = m3$param[4]
m3_scale_lat = m3$param[5]
m3_scale_lon = m3$param[6]

m3_shape_intercept = m3$param[7]
m3_shape_lat = m3$param[8]
m3_shape_lon = m3$param[9]

### d) Modèle 4 --------
loc.form4 <- ~ lat + lon
scale.form4 <- ~ lat 
shape.form4 <- ~ lat + lon

m4 = fitspatgev(D_Mat, coord, loc.form4, scale.form4, shape.form4, corr =TRUE)
#m2 <- fitspatgev(t(PL[,-c(1,2)]), coord, loc.form1, scale.form1, shape.form1)

m4_loc_intercept = m4$param[1]
m4_loc_lat = m4$param[2]
m4_loc_lon = m4$param[3]

m4_scale_intercept = m4$param[4]
m4_scale_lat = m4$param[5]

m4_shape_intercept = m4$param[6]
m3_shape_lat = m4$param[7]
m3_shape_lon = m4$param[8]

## e) Modèle Mod3 ---------
years <- rep(1:24, each = ncol(D_Mat))
temp.cov = matrix(years, 24, byrow=T)
colnames(temp.cov) =c(1:1568)

loc.form_3 <- y ~ lat + lon
scale.form_3 <- y ~ lat + lon
shape.form_3 <- y ~ 1
temp.form.loc_3 <- ~ poly(years, 1)
temp.form.scale_3 <- ~ poly(years, 1)

Mod3 <- fitspatgev(D_Mat, coord, loc.form_3, scale.form_3, shape.form_3, 
                   temp.form.loc = temp.form.loc_3, temp.form.scale = temp.form.scale_3, 
                   temp.cov = temp.cov, corr =TRUE)
    
Mod3_loc_intercept = Mod3$param[1]
Mod3_loc_lat = Mod3$param[2]
Mod3_loc_lon = Mod3$param[3]
Mod3_scale_intercept = Mod3$param[4]
Mod3_scale_lat = Mod3$param[5]
Mod3_scale_lon = Mod3$param[6]
Mod3_shape_intercept = Mod3$param[7]
Mod3_temp_loc = Mod3$param[8]
Mod3_temp_scale = Mod3$param[9]


## Critères ----------

Deviance_m1 = m1$deviance
Deviance_m2 = m2$deviance
Deviance_m3 = m3$deviance
Deviance_m4 = m4$deviance

# BIC
bic_funct <- function(modele, n) {
  k <- length(modele$param)
  bic_mod <- -2 * modele$logLik + k * log(n)  # n = nombre d'observations
  return(bic_mod)
}

n = 24 * ncol(D_Mat)
BIC1 = bic_funct(m1, n)
BIC2 = bic_funct(m2, n)
BIC3 = bic_funct(m3, n)
BIC4 = bic_funct(m4, n)

# AIC
aic_funct = function(modele){
  k = length(modele$param)
  aic_mod = -2*modele$logLik + 2*k
  cat("\nlogLik =", modele$logLik, "; k = ", k)
  cat("\nAIC = ", aic_mod)
  return(aic_mod)
}
AIC1 = aic_funct(m1)
AIC2 = aic_funct(m2)
AIC3 = aic_funct(m3)
AIC4 = aic_funct(m4)

# Anova
anova(m1, m3)

# Tableau des critères
tab_critere = data.frame(#Modeles = c("m1", "m2", "m3", "m4"; "Mod3"),
                         Deviances = c(Deviance_m1, Deviance_m2, Deviance_m3, Deviance_m4, Mod3$deviance),
                         AIC = c(AIC1, AIC2, AIC3, AIC4, aic_funct(Mod3)),
                         BIC = c(BIC1, BIC2, BIC3, BIC4, bic_funct(Mod3, n)),
                         TIC = c(TIC(m1), TIC(m2), TIC(m3), TIC(m4), TIC(Mod3))
                         )
tab_critere
write.csv2(tab_critere, 'tab_critère_classic.csv')

# 2. Modèles max-stables -----------

D <- dist(coord)
weights <- 1 / (1 + D)

start.trend = as.list(m1$fitted.values)

loc.form_3 <- y ~ lat + lon
scale.form_3 <- y ~ lat + lon
shape.form_3 <- y ~ 1
#temp.form.loc_3 <- ~ poly(years, 1)
#temp.form.scale_3 <- ~ poly(years, 1)
##a0) Modèle de smith ------------
start = c(list(cov11 = var(coord[,1]), cov12 = cov(coord[,1],coord[,2]), cov22 = var(coord[,2])), start.trend)
smith <- fitmaxstab(D_Mat, coord, "gauss", loc.form_3, scale.form_3, shape.form_3, 
                     #temp.cov = temp.cov, temp.form.loc = temp.form.loc_3, temp.form.scale = temp.form.scale_3,
                     method = "nlm", weights = weights,
                     start = start, typsize = unlist(start))

## a) Modèle de Schlater --------
start = c(list(range = 5, smooth = 1.345), start.trend)
schlat <- fitmaxstab(D_Mat, coord, "powexp", loc.form_3, scale.form_3, shape.form_3, 
                     #temp.cov = temp.cov, temp.form.loc = temp.form.loc_3, temp.form.scale = temp.form.scale_3,
                     nugget = 0, method = "nlm", weights = weights,
                     start = start, typsize = unlist(start))

schlat_nugget = schlat$param[1] # 

schlat_range = schlat$param[2] # lambda

schlat_smooth =schlat$param[3] # Kappa

schlat_loc_intercept = schlat$param[4]
schlat_loc_lat = schlat$param[5]
schlat_loc_lon = schlat$param[6]

schlat_scale_intercept = schlat$param[7]
schlat_scale_lat = schlat$param[8]
schlat_scale_lon = schlat$param[9]

schlat_shape_intercept = schlat$param[10]

schlat_loc_temp = schlat$param[11]
schlat_scale_temp = schlat$param[12]

schlat_dev = schlat$deviance
schlat_opt.val = schlat$opt.value


schlat_estimates = as.data.frame(schlat$fitted.values)
write.csv2(schlat_estimates, 'schlat_estimate_t.csv')


# h = 
extcoeff(schlat, 'powexp', 200)

# Compare the empirical F-madogram cloud to the fitted extremal coefficient function
concprob(D_Mat, coord, schlat, n.bins = 100, col = c("black", "red"), which = "emp")
concprob(D_Mat, coord, schlat, n.bins = 500, which = "boot")
concprob(D_Mat, coord, schlat, n.bins = 10000, col = c("gray", "red"), which = "kendall", compute.std.err = TRUE, ylim = c(0, 2))


## b) Modèle de Brown-Resnick --------
start <- c(list(range = 1.71, smooth = 0.46), start.trend)

brown_t <- fitmaxstab(D_Mat, coord, "brown", loc.form_3, scale.form_3, shape.form_3, 
                      #temp.cov = temp.cov, temp.form.loc = temp.form.loc_3, temp.form.scale = temp.form.scale_3,
           method = "nlm", weights = weights,
           start = start, typsize = unlist(start))

Brown_estimates = as.data.frame(brown$fitted.values)
write.csv2(Brown_estimates, 'Brown_estimate_t.csv')


brown_range = brown_t$param[1]

brown_smooth =brown_t$param[2]

brown_loc_intercept = brown_t$param[3]
brown_loc_lat = brown_t$param[4]
brown_loc_lon = brown_t$param[5]

brown_scale_intercept = brown$param[6]
brown_scale_lat = brown$param[7]
brown_scale_lon = brown$param[8]

brown_shape_intercept = brown$param[9]

brown_dev = brown$deviance
brown_opt.val = brown$opt.value

# Matrice hessienne
brown_ihessian = as.data.frame(brown$ihessian)
write.csv2(brown_ihessian, 'brown_matrice_hessienne.csv')

# Score de la variance
brown_var.score = as.data.frame(brown$var.score)
write.csv2(brown_var.score, 'brown_matrice_var.score.csv')

# Variance covariance
brown_var.cov = as.data.frame(brown$var.cov)
write.csv2(brown_var.cov, 'brown_matrice_var.cov.csv')

brown_estimates = as.data.frame(brown$fitted.values)
brown_std.error = as.data.frame(brown$std.err)
brown_estimates$brown_std.error = brown_std.error$`brown$std.err`
write.csv2(brown_estimates, 'brown_estimate_std.error.csv')
# h = 
extcoeff(brown, 'brown', c(brown_nugget, brown_range, brown_smooth), 200)

# Compare the empirical F-madogram cloud to the fitted extremal coefficient function
concprob(D_Mat, coord, brown, n.bins = 500, col = c("black", "red"), which = "emp", add = TRUE)
concprob(D_Mat, coord, brown, n.bins = 500, which = "boot")
concprob(D_Mat, coord, brown, n.bins = 10000, col = c("gray", "red"), which = "kendall", compute.std.err = TRUE, xlim = c(0, 40), ylim = c(0, 1))




## c) Modèle Extremal-t --------
start <- c(list(range = 100, smooth = 0.53, DoF = 7.45), start.trend)
extt_t <- fitmaxstab(D_Mat, coord, "tpowexp", loc.form_3, scale.form_3, shape.form_3, 
                     #temp.cov = temp.cov, temp.form.loc = temp.form.loc_3, temp.form.scale = temp.form.scale_3,
           nugget = 0, method = "nlm", weights = weights,
           start = start, typsize = unlist(start))

extt_t_estimates = as.data.frame(extt_t$fitted.values)
write.csv2(extt_t_estimates, 'extt_estimate_t.csv')

### c.0. extt1-------------
fmad0 <- fmadogram(D_Mat, coord, n.bins = 10000) 

# construire poids binaire basé sur distance
D <- dist(coord)
#weights <- ifelse(D <= 4, 1, 0)
weights <- 1 / (1 + D)


start2 <- c(list(range = 2, smooth = 0.5, DoF = 7.5), start.trend)
extt2 <- fitmaxstab(D_Mat, coord, "tpowexp", loc.form_3, scale.form_3, shape.form_3, 
                    #temp.cov = temp.cov, temp.form.loc = temp.form.loc_3, temp.form.scale = temp.form.scale_3,
                    nugget = 0, method = "nlm", weights = weights,
                    start = start, typsize = unlist(start))

extt1_nugget = extt1$param[1]

extt1_range = extt1$param[2]

extt1_smooth =extt1$param[3]

extt1_DoF =extt1$param[4]

extt1_loc_intercept = extt1$param[5]
extt1_loc_lat = extt1$param[6]
extt1_loc_lon = extt1$param[7]

extt1_scale_intercept = extt1$param[8]
extt1_scale_lat = extt1$param[9]
extt1_scale_lon = extt1$param[10]

extt1_shape_intercept = extt1$param[11]

extt1_dev = extt1$deviance
extt1_opt.val = extt1$opt.value
extt1_tic = TIC(extt1)



# Fun cov
extt1_cov.fun = as.data.frame(extt1$cov.fun)
write.csv2(extt1_cov.fun, 'extt1_cov.fun.csv')

# Variance cov
extt1_var.cov = as.data.frame(extt1$var.cov)
write.csv2(extt1_var.cov, 'extt1_matrice_var.cov.csv')

# Graphe de la matrice des corrélation
mat_cov1 <- as.matrix(extt1_var.cov)
matrice_cor1 <- cov2cor(mat_cov1) # transformation en matrice de corr
corrplot(
  matrice_cor1, 
  method = "shade",      # Type de représentation : "circle", "square", "number", "shade"
  type = "upper",         # Afficher seulement la partie supérieure (évite la redondance)
  order = "hclust",       # Réorganiser les variables par clustering hiérarchique
  tl.col = "black",       # Couleur des étiquettes (noms des variables)
  tl.srt = 45,            # Rotation des étiquettes
  diag = TRUE            # Ne pas afficher les coefficients sur la diagonale
)

# Matrice hessienne
extt1_ihessian = as.data.frame(extt1$ihessian)
write.csv2(extt1_ihessian, 'extt1_matrice_hessienne.csv')

# Score de la variance
extt1_var.score = as.data.frame(extt1$var.score)
write.csv2(extt1_var.score, 'extt1_matrice_var.score.csv')

extt1_estimates = as.data.frame(extt1$fitted.values)
extt1_std.error = as.data.frame(extt1$std.err)
extt1_estimates$extt1_std.error = extt1_std.error$`extt1$std.err`
write.csv2(extt1_estimates, 'extt1_estimate_std.error.csv')



fmad1 <- fmadogram(D_Mat, coord, extt1, col = c("gray", "red"), n.bins = 10000)  
write.csv2(fmad1, "fmadogramm1bins10000.csv")
#plot(fmad1)  # montre F-madogram et binned estimates
legend("top", c("Empirique","Modèle"), col=c("gray","red"))
# convertir en extremal coefficient via relation: theta = (1+2*v)/(1-2*v)
#v <- fmad1[,2] # or appropriate element depending on returned structure
#theta_emp1 <- (1 + 2*v)/(1 - 2*v)
extt_concprob1 = concprob(D_Mat, coord, extt1, col = c("gray", "red"), which = "kendall", compute.std.err = TRUE, ylim = c(0, 1))
write.csv2(extt_concprob1, "extt1_concprob.csv")


### c.1. Extraction des paramètres---------
extt_nugget = extt$param[1]

extt_range = extt$param[2]

extt_smooth =extt$param[3]

extt_DoF =extt$param[4]

extt_loc_intercept = extt$param[5]
extt_loc_lat = extt$param[6]
extt_loc_lon = extt$param[7]

extt_scale_intercept = extt$param[8]
extt_scale_lat = extt$param[9]
extt_scale_lon = extt$param[10]

extt_shape_intercept = extt$param[11]

extt_dev = extt$deviance
extt_opt.val = extt$opt.value
extt_tic = TIC(extt)

# Function cov
# Fun cov
extt_cov.fun = as.data.frame(extt$cov.fun)
write.csv2(extt_cov.fun, 'extt_cov.fun.csv')

# Variance cov
extt_var.cov = as.data.frame(extt$var.cov)
write.csv2(extt_var.cov, 'extt_matrice_var.cov.csv')

# Graphe de la matrice des corrélation
mat_cov <- as.matrix(extt_var.cov)
matrice_cor <- cov2cor(mat_cov) # transformation en matrice de corr
corrplot(
  matrice_cor, 
  method = "shade",      # Type de représentation : "circle", "square", "number", "shade"
  type = "upper",         # Afficher seulement la partie supérieure (évite la redondance)
  order = "hclust",       # Réorganiser les variables par clustering hiérarchique
  tl.col = "black",       # Couleur des étiquettes (noms des variables)
  tl.srt = 45,            # Rotation des étiquettes
  diag = TRUE            # Ne pas afficher les coefficients sur la diagonale
)

# Matrice hessienne
extt_ihessian = as.data.frame(extt$ihessian)
write.csv2(extt_ihessian, 'extt_matrice_hessienne.csv')

# Score de la variance
extt_var.score = as.data.frame(extt$var.score)
write.csv2(extt_var.score, 'extt_matrice_var.score.csv')

extt_estimates = as.data.frame(extt$fitted.values)
extt_std.error = as.data.frame(extt$std.err)
extt_estimates$extt_std.error = extt_std.error$`extt$std.err`
write.csv2(extt_estimates, 'extt_estimate_std.error.csv')



### c.2. Diagnostic et validation -----------
# Critère de choix du meilleur modèle
tab_critere_max_stable = data.frame(Modeles = c("Schlater", "BR", "Extremal-t", "Extremal-t1"),
                                    Deviances = c(schlat$deviance, brown$deviance, extt$deviance, extt1$deviance),
                                    TIC = c("NA", TIC(brown), TIC(extt), TIC(extt1)),
                                    logLik = c(schlat$logLik, brown$logLik, extt$logLik, extt1$logLik),
                                    Opt.value = c(schlat$opt.value, brown$opt.value, extt$opt.value, extt1$opt.value)
)
tab_critere_max_stable
write.csv2(tab_critere_max_stable, 'maxstable_critères.csv')

tab_critere_max_stable_t = data.frame(Deviances = c(schlat$deviance, brown_t$deviance, extt_t$deviance),
                                      TIC = c(TIC(schlat), TIC(brown_t), TIC(extt_t)),
                                      logLik = c(schlat$logLik, brown_t$logLik, extt_t$logLik),
                                      Opt.value = c(schlat$opt.value, brown_t$opt.value, extt_t$opt.value)
)
tab_critere_max_stable_t
write.csv2(tab_critere_max_stable_t, 'maxstable_critères_t.csv')

# MEILLEUR MODEL EXTREMAL-T


# fmadogram et coef extremal en fonction de la distance en dégré
fmad <- fmadogram(D_Mat, coord, extt_t, col = c("gray", "red"), n.bins = 10000)  
write.csv2(fmad, "fmadogramm1000bins.csv")
#plot(fmad)  # montre F-madogram et binned estimates
legend("topright", c("Empirique","Modèle"), col=c("gray","red"))
# convertir en extremal coefficient via relation: theta = (1+2*v)/(1-2*v)
#v <- fmad[,2] # or appropriate element depending on returned structure
#theta_emp <- (1 + 2*v)/(1 - 2*v)

# prédire theta(h) à partir du modèle 
#h_D = as.integer(max(D))
# SpatialExtremes propose extcoeff() ou fonctions équivalentes
#theta_fit <- extcoeff(extt, distance = seq(0, h_D, length = 50))
# tracé comparatif

#plot(NULL, xlim = c(0, h_D), ylim = c(1,2), xlab="h", ylab=TeX("$theta (h)$"))
#points(fmad[,1], theta_emp, pch=19, col="gray")
#lines(tab_coef_extremal$h, tab_coef_extremal$extt_theta.coef, col="red")
#lines(theta_fit$dist, theta_fit$theta, col="red", lwd=2)
#legend("topright", c("Empirique","Modèle"), col=c("gray","red"), pch=c(19,NA), lty=c(NA,1))



# Cette fonction calcule les estimations empiriques ou extrêmes de la probabilité de concurrence par paires.
concprob(D_Mat, coord, extt_t, n.bins = 500, col = c("black", "red"), which = "emp", add = TRUE)
concprob(D_Mat, coord, extt_t, n.bins = 500, which = "boot")
extt_concprob = concprob(D_Mat, coord, extt_t, col = c("gray", "red"),
                         which = "kendall", compute.std.err = TRUE, ylim = c(0, 1),
                         plot = TRUE)
plot(extt_concprob[,1], extt_concprob[,2])
write.csv2(extt_concprob, "extt_concprob.csv")

h_D=seq(0, 5, 0.5)
# Theta de h en fonction de h (1: dépendance parfaite; 2: indépendance)
tab_coef_extremal = data.frame(
  h = h_D,
  schlat_theta.coef = schlat$ext.coeff(h_D),
  brown_theta.coef = brown_t$ext.coeff(h_D),
  extt_theta.coef = extt_t$ext.coeff(h_D),
  extt_cov.fun = extt$cov.fun(h_D) # function covariance
)
tab_coef_extremal
write.csv2(tab_coef_extremal, 'coef_extremal_theta_h.csv')

plot(tab_coef_extremal$h, tab_coef_extremal$extt_theta.coef, xlab = "h", ylab = TeX("$theta (h)$"))
lines(tab_coef_extremal$h, rep(1.5, length(h_D)))
plot(tab_coef_extremal$h, tab_coef_extremal$brown_theta.coef, xlab = "h", ylab = TeX("$theta (h)$"))
plot(tab_coef_extremal$h, tab_coef_extremal$schlat_theta.coef, xlab = "h", ylab = TeX("$theta (h)$"))
plot(tab_coef_extremal$h, tab_coef_extremal$extt_cov.fun, xlab = "h", ylab = TeX("$rho (h)$"))


# 3. Diagnostic complet --------
## 3.0 DIAGNOSTIC COMPLET DU MODÈLE ---------
#
library(SpatialExtremes)
library(ggplot2)
library(dplyr)
library(gridExtra)

## 3.1. INFORMATIONS GÉNÉRALES DU MODÈLE-----------

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║          DIAGNOSTIC DU MODÈLE MAX-STABLE                      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Résumé du modèle
cat("=== INFORMATIONS DU MODÈLE ===\n")
cat("Type de modèle: T-Power Exponential\n")
cat("Nombre de sites:", ncol(D_Mat), "\n")
cat("Nombre d'années:", nrow(D_Mat), "\n")
cat("Nombre total d'observations:", nrow(D_Mat) * ncol(D_Mat), "\n")
cat("Méthode d'estimation: Pairwise likelihood (weights)\n")
cat("Convergence:", extt_t$convergence, "(0 = OK)\n\n")

# Paramètres estimés
cat("=== PARAMÈTRES ESTIMÉS ===\n")
print(round(extt_t$param, 4))
cat("\n")

# Erreurs standard (si disponibles)
if(!is.null(extt_t$std.err)) {
  cat("=== ERREURS STANDARD ===\n")
  print(round(extt_t$std.err, 4))
  cat("\n")
}

## 3.2. CALCUL DU CLIC----------

calcul_clic <- function(model) {
  # Log-vraisemblance composite (négative car optimisation minimise)
  loglik_composite <- -model$opt.value
  
  # Nombre de paramètres
  n_params <- length(model$param)
  
  # Nombre de sites
  n_sites <- ncol(model$data)
  
  # CLIC
  clic <- -2 * loglik_composite + 2 * n_params
  
  # CLIC corrigé pour petits échantillons
  if(n_sites > n_params + 1) {
    clicc <- clic + (2 * n_params * (n_params + 1)) / (n_sites - n_params - 1)
  } else {
    clicc <- NA
  }
  
  # BIC version composite
  bic_composite <- -2 * loglik_composite + n_params * log(n_sites)
  
  return(list(
    loglik = loglik_composite,
    CLIC = clic,
    CLICc = clicc,
    BIC = bic_composite,
    n_params = n_params,
    n_sites = n_sites
  ))
}

calcul_clic(smith)

clic_results <- calcul_clic(extt_t)

cat("=== CRITÈRES D'INFORMATION ===\n")
cat("Log-vraisemblance composite:", round(clic_results$loglik, 2), "\n")
cat("CLIC:", round(clic_results$CLIC, 2), "\n")
cat("CLICc:", round(clic_results$CLICc, 2), "\n")
cat("BIC:", round(clic_results$BIC, 2), "\n")
cat("Nombre de paramètres:", clic_results$n_params, "\n\n")

## 3.3. CALCUL DES MESURES EMPIRIQUES------

cat("=== CALCUL DES MESURES EMPIRIQUES ===\n")
cat("Calcul du coefficient d'extrémité empirique...\n")
emp_ext <- fmadogram(D_Mat, coord, which = "ext")

cat("Calcul du F-madogram empirique...\n")
emp_F <- fmadogram(D_Mat, coord, which = "F")

cat("Calcul du variogramme empirique...\n")
emp_vario <- fmadogram(D_Mat, coord, which = "vario")

cat("Calculs terminés.\n\n")

## 3.4. EXTRACTION DU COEFFICIENT D'EXTRÉMITÉ THÉORIQUE------

cat("=== EXTRACTION DU COEFFICIENT THÉORIQUE ===\n")

# Créer une séquence de distances
max_dist <- max(emp_ext[,1])
distances_seq <- seq(0.01, max_dist, length.out = 500)

# Fonction pour calculer le coefficient d'extrémité théorique
# Pour le modèle t-power exponential
calcul_ext_coef_theorique <- function(h, range, smooth, DoF) {
  # h: distance
  # range: paramètre de portée
  # smooth: paramètre de lissage
  # DoF: degrés de liberté
  
  if(h == 0) return(1)
  
  # Coefficient d'extrémité pour t-power exponential
  # θ(h) = 2 * Φ(sqrt((h/range)^smooth / 2))
  # où Φ est la fonction de répartition de Student-t
  
  arg <- sqrt((h / range)^smooth / 2)
  theta <- 2 * pt(arg * sqrt(DoF), df = DoF)
  
  return(theta)
}

# Pour le modèle brown 
calcul_brown_coef_theorique <- function(h, range, smooth) {
  # h: distance
  # range: paramètre de portée
  # smooth: paramètre de lissage
  
  if(h == 0) return(1)
  #Semi-varogramme
  gamma_h = (h/range)^smooth
  
  # Coefficient d'extrémité pour brown
  # θ(h) = 2 * Φ(sqrt(gamma(h) / 2))
  # où Φ est la fonction de répartition de la loi normale
  
  
  theta <- 2 * pnorm(sqrt(gamma_h/2))
  
  return(theta)
}

calcul_schlat_coef_theorique <- function(h, range, smooth) {
  # h: distance
  # range: paramètre de portée
  # smooth: paramètre de lissage
  
  if(h == 0) return(1)
  #Semi-varogramme
  rho_h = exp(-(h/range)^smooth)
  
  # Coefficient d'extrémité pour brown
  
  
  theta <- 1 + 1/2*sqrt(1+2/rho_h)
  
  return(theta)
}

calcul_smith_coef_theorique <- function(h, sigma) {
  # h: distance
  # range: paramètre de portée
  # smooth: paramètre de lissage
  
  if(h == 0) return(1)
  sigma_inv = solve(sigma) 
  
  mahal_sq = dist^2 * sigma_inv[1,1]
  mahal = sqrt(mahal_sq)
  theta <- 2*pnorm(mahal/2)
  
  return(theta)
}


# Extraire les paramètres du modèle
range_param <- extt_t$param["range"]
smooth_param <- extt_t$param["smooth"]
DoF_param <- extt_t$param["DoF"]

cat("Paramètres de dépendance spatiale:\n")
cat("  Range:", round(range_param, 4), "\n")
cat("  Smooth:", round(smooth_param, 4), "\n")
cat("  DoF:", round(DoF_param, 4), "\n\n")

# Calculer le coefficient théorique pour toutes les distances
ext_coef_theorique <- sapply(distances_seq, function(h) {
  calcul_ext_coef_theorique(h, range_param, smooth_param, DoF_param)
})

brown_ext_coef_theorique <- sapply(distances_seq, function(h) {
  calcul_brown_coef_theorique(h, brown_t$param["range"], brown_t$param["smooth"])
})

schlat_ext_coef_theorique <- sapply(distances_seq, function(h) {
  calcul_schlat_coef_theorique(h, schlat$param["range"], schlat$param["smooth"])
})

smith_ext_coef_theorique <- sapply(distances_seq, function(h) {
  calcul_smith_coef_theorique(h, smith$var.cov)
})

# Créer un data.frame
df_theorique <- data.frame(
  distance = distances_seq,
  ext_coef = ext_coef_theorique,
  brown_ext_coef = brown_ext_coef_theorique,
  schlat_ext_coef = schlat_ext_coef_theorique
)



## 3.5. IDENTIFIER LA DISTANCE h POUR θ(h) ∈ [1.3, 1.7]-----

cat("=== IDENTIFICATION DE LA DISTANCE h ===\n")
cat("Recherche de h tel que θ(h) ∈ [1.3, 1.7]\n\n")

# Trouver les distances correspondantes
distances_cible <- df_theorique %>%
  filter(ext_coef >= 1.3 & ext_coef <= 1.7)

if(nrow(distances_cible) > 0) {
  h_min <- min(distances_cible$distance)
  h_max <- max(distances_cible$distance)
  h_median <- median(distances_cible$distance)
  
  theta_at_h_min <- distances_cible$ext_coef[which.min(distances_cible$distance)]
  theta_at_h_max <- distances_cible$ext_coef[which.max(distances_cible$distance)]
  theta_at_h_median <- distances_cible$ext_coef[which.min(abs(distances_cible$distance - h_median))]
  
  cat("╔════════════════════════════════════════════════════════════════╗\n")
  cat("║  RÉSULTATS: Distance h pour θ(h) ∈ [1.3, 1.7]                ║\n")
  cat("╚════════════════════════════════════════════════════════════════╝\n")
  cat(sprintf("Distance minimale (h_min):  %.4f  →  θ(h) = %.4f\n", h_min, theta_at_h_min))
  cat(sprintf("Distance médiane (h_med):   %.4f  →  θ(h) = %.4f\n", h_median, theta_at_h_median))
  cat(sprintf("Distance maximale (h_max):  %.4f  →  θ(h) = %.4f\n", h_max, theta_at_h_max))
  cat(sprintf("\nIntervalle: [%.4f, %.4f]\n", h_min, h_max))
  cat(sprintf("Largeur de l'intervalle: %.4f\n\n", h_max - h_min))
  
  # Valeurs spécifiques
  cat("Valeurs spécifiques:\n")
  
  # θ(h) = 1.3
  idx_1.3 <- which.min(abs(df_theorique$ext_coef - 1.3))
  h_1.3 <- df_theorique$distance[idx_1.3]
  cat(sprintf("  Pour θ(h) = 1.3  →  h ≈ %.4f\n", h_1.3))
  
  # θ(h) = 1.5
  idx_1.5 <- which.min(abs(df_theorique$ext_coef - 1.5))
  h_1.5 <- df_theorique$distance[idx_1.5]
  cat(sprintf("  Pour θ(h) = 1.5  →  h ≈ %.4f\n", h_1.5))
  
  # θ(h) = 1.7
  idx_1.7 <- which.min(abs(df_theorique$ext_coef - 1.7))
  h_1.7 <- df_theorique$distance[idx_1.7]
  cat(sprintf("  Pour θ(h) = 1.7  →  h ≈ %.4f\n\n", h_1.7))
  
} else {
  cat("ATTENTION: Aucune distance ne correspond à θ(h) ∈ [1.3, 1.7]\n")
  cat("Plage actuelle de θ(h): [", round(min(ext_coef_theorique), 2), ", ", 
      round(max(ext_coef_theorique), 2), "]\n\n")
}

## 3.6. CALCUL DES MÉTRIQUES DE VALIDATION-------

cat("=== MÉTRIQUES DE VALIDATION ===\n")

# Interpoler le modèle théorique aux distances empiriques
ext_coef_model_at_emp <- sapply(emp_ext[,1], function(h) {
  calcul_ext_coef_theorique(h, range_param, smooth_param, DoF_param)
})

ext_coef_model_at_emp_brown <- sapply(emp_ext[,1], function(h) {
  calcul_brown_coef_theorique(h, brown_t$param["range"], brown_t$param["smooth"])
})

ext_coef_model_at_emp_schlat <- sapply(emp_ext[,1], function(h) {
  calcul_schlat_coef_theorique(h, schlat$param["range"], schlat$param["smooth"])
})

ext_coef_model_at_emp_smith <- sapply(emp_ext[,1], function(h) {
  calcul_smith_coef_theorique(h, smith$var.cov)
})


# Calculer les métriques
rmse_ext <- sqrt(mean((emp_ext[,3] - ext_coef_model_at_emp)^2, na.rm = TRUE))
mae_ext <- mean(abs(emp_ext[,3] - ext_coef_model_at_emp), na.rm = TRUE)
cor_ext <- cor(emp_ext[,3], ext_coef_model_at_emp, use = "complete.obs")
biais_ext = mean(ext_coef_model_at_emp - emp_ext[,3], na.rm = TRUE)

rmse_ext_br <- sqrt(mean((emp_ext[,3] - ext_coef_model_at_emp_brown)^2, na.rm = TRUE))
mae_ext_br <- mean(abs(emp_ext[,3] - ext_coef_model_at_emp_brown), na.rm = TRUE)
cor_ext_br <- cor(emp_ext[,3], ext_coef_model_at_emp_brown, use = "complete.obs")
biais_ext_br = mean(ext_coef_model_at_emp_brown - emp_ext[,3], na.rm = TRUE)

rmse_ext_schl <- sqrt(mean((emp_ext[,3] - ext_coef_model_at_emp_schlat)^2, na.rm = TRUE))
mae_ext_schl <- mean(abs(emp_ext[,3] - ext_coef_model_at_emp_schlat), na.rm = TRUE)
cor_ext_schl <- cor(emp_ext[,3], ext_coef_model_at_emp_schlat, use = "complete.obs")
biais_ext_schl = mean(ext_coef_model_at_emp_schlat - emp_ext[,3], na.rm = TRUE)

rmse_ext_smith <- sqrt(mean((emp_ext[,3] - ext_coef_model_at_emp_smith)^2, na.rm = TRUE))
mae_ext_smith <- mean(abs(emp_ext[,3] - ext_coef_model_at_emp_smith), na.rm = TRUE)
cor_ext_smith <- cor(emp_ext[,3], ext_coef_model_at_emp_smith, use = "complete.obs")
biais_ext_smith = mean(ext_coef_model_at_emp_smith - emp_ext[,3], na.rm = TRUE)


cat("Coefficient d'extrémité:\n")
cat(sprintf("  RMSE: %.4f\n", rmse_ext))
cat(sprintf("  MAE:  %.4f\n", mae_ext))
cat(sprintf("  COR:  %.4f\n\n", cor_ext))

## 3.7. VISUALISATIONS---------

cat("=== GÉNÉRATION DES GRAPHIQUES ===\n")

# Graphique 1: Coefficient d'extrémité avec zone d'intérêt
p1 <- ggplot() +
  # Points empiriques
  geom_point(data = data.frame(dist = emp_ext[,1], 
                               ext = emp_ext[,3]),
             aes(x = dist, y = ext), 
             color = "gray50", size = 2.5, alpha = 0.6) +
  # Courbe théorique
  geom_line(data = df_theorique, 
            aes(x = distance, y = ext_coef), 
            color = "red", size = 1.2) +
  # Zone d'intérêt [1.3, 1.7]
  geom_hline(yintercept = c(1.3, 1.7), 
             linetype = "dashed", color = "darkgreen", size = 0.8) +
  annotate("rect", xmin = 0, xmax = max_dist, 
           ymin = 1.3, ymax = 1.7,
           alpha = 0.1, fill = "green") +
  # Limites théoriques
  geom_hline(yintercept = 1, linetype = "dotted", color = "black") +
  geom_hline(yintercept = 2, linetype = "dotted", color = "black") +
  # Labels
  annotate("text", x = max_dist * 0.9, y = 1, 
           label = "Perfect dependence", size = 3, color = "black") +
  annotate("text", x = max_dist * 0.9, y = 2, 
           label = "Independence", size = 3, color = "black") +
  annotate("text", x = max_dist * 0.5, y = 1.5, 
           label = "\n[1.3, 1.7]", 
           size = 4, color = "darkgreen", fontface = "bold") +
  labs(title = "Extremity Coefficient: Model vs Empirical",
       subtitle = paste0("RMSE = ", round(rmse_ext, 3), 
                         " | BIAS = ", round(biais_ext, 3),
                         " | Correlation = ", round(cor_ext, 3)),
       x = "Distance h (in degree)", 
       y = "Extremal Coefficient θ(h)",
       caption = "Red = Theoretical model | Gray = Empirical") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank())

# Si on a trouvé les distances
if(exists("h_1.3") && exists("h_1.5") && exists("h_1.7")) {
  p1 <- p1 +
    geom_vline(xintercept = c(h_1.3, h_1.5, h_1.7), 
               linetype = "dashed", color = "blue", alpha = 0.5) +
    annotate("text", x = h_1.3, y = 2.1, 
             label = sprintf("h(1.3)=%.2f", h_1.3), 
             angle = 90, size = 3, color = "blue") +
    annotate("text", x = h_1.5, y = 2.1, 
             label = sprintf("h(1.5)=%.2f", h_1.5), 
             angle = 90, size = 3, color = "blue") +
    annotate("text", x = h_1.7, y = 2.1, 
             label = sprintf("h(1.7)=%.2f", h_1.7), 
             angle = 90, size = 3, color = "blue")
}

print(p1)

# Graphique 2: Zoom sur la zone d'intérêt
if(nrow(distances_cible) > 0) {
  p1_zoom <- ggplot() +
    
    geom_point(data = data.frame(dist = emp_ext[,1], 
                                 ext = emp_ext[,3]) %>%
                 filter(ext >= 1.2 & ext <= 1.8),
               aes(x = dist, y = ext), 
               color = "gray50", size = 3, alpha = 0.7) +
    
    geom_line(data = df_theorique %>% filter(ext_coef >= 1.2 & ext_coef <= 1.8), 
              aes(x = distance, y = ext_coef), 
              color = "red", size = 1.5) +
    
    geom_hline(yintercept = c(1.3, 1.7), 
               linetype = "dashed", color = "darkgreen", size = 1) +
    annotate("rect", xmin = h_min, xmax = h_max, 
             ymin = 1.25, ymax = 1.75,
             alpha = 0.15, fill = "green") +
    labs(title = "Zoom: Zone θ(h) ∈ [1.3, 1.7]",
         subtitle = sprintf("Distance interval: [%.3f, %.3f]", h_min, h_max),
         x = "Distance", 
         y = "θ(h)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  print(p1_zoom)
}

# Graphique 3: F-madogram
p2 <- ggplot() +
  # geom_line(data = data.frame(x = seq(0, max(emp_F[,1]), length.out = 100)),
  #       aes(x = x), color = "blue", size = 1, stat = "identity") +
  geom_point(data = data.frame(dist = emp_F[,1], F = emp_F[,2]),
             aes(x = dist, y = F), 
             color = "red", size = 2, alpha = 0.6) +
  labs(title = "F-madogram: Model vs Empirical",
       x = "Distance", y = "F-madogram") +
  theme_minimal(base_size = 12)

# Utiliser la fonction plot de SpatialExtremes pour le F-madogram
plot(extt, which = "F", main = "F-madogram: Model vs Empirical",
     col = "red", lwd = 2, xlab = "Distance", ylab = "F")
points(emp_F[,1], emp_F$F, pch = 16, col = "gray50", cex = 1.2)
legend("topright", c("Model", "Empirical"), 
       col = c("red", "gray50"), lty = c(1, NA), pch = c(NA, 16), lwd = 2)

# Graphique 4: Variogramme
plot(extt, which = "vario", main = "Variogram: Model vs Empirical",
     col = "red", lwd = 2, xlab = "Distance", ylab = "Semi-variance")
points(emp_vario[,1], emp_vario$vario, pch = 16, col = "gray50", cex = 1.2)
legend("bottomright", c("Model", "Empirical"), 
       col = c("red", "gray50"), lty = c(1, NA), pch = c(NA, 16), lwd = 2)

## 3.8. PANNEAU DE DIAGNOSTICS COMPLET--------

pdf("diagnostic_maxstab_complet.pdf", width = 16, height = 10)

par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 2))

# 1. Coefficient d'extrémité
plot(df_theorique$distance, df_theorique$ext_coef, 
     type = "l", col = "red", lwd = 2.5,
     main = "Extremal coefficient",
     xlab = "Distance", ylab = "theta(h)",
     ylim = c(1, 2))
points(emp_ext[,1], emp_ext[,3], pch = 16, col = "gray50", cex = 1.2)
abline(h = c(1, 1.3, 1.7, 2), lty = c(3, 2, 2, 3), col = c("black", "darkgreen", "darkgreen", "black"))
if(exists("h_min")) {
  abline(v = c(h_min, h_max), lty = 2, col = "blue")
}
legend("bottomright", c("Model", "Empirical", "Zone [1.3, 1.7]"), 
       col = c("red", "gray50", "darkgreen"), 
       lty = c(1, NA, 2), pch = c(NA, 16, NA), lwd = c(2, NA, 1),
       cex = 0.8, bty = "n")

# 2. F-madogram
plot(extt, which = "F", col = "blue", lwd = 2.5,
     main = "F-madogram", xlab = "Distance", ylab = "F")
points(emp_F$dist, emp_F$F, pch = 16, col = "red", cex = 1.2)
legend("topright", c("Model", "Empirical"), 
       col = c("blue", "red"), lty = c(1, NA), pch = c(NA, 16), 
       lwd = 2, cex = 0.8, bty = "n")

# 3. Variogramme
plot(extt, which = "vario", col = "blue", lwd = 2.5,
     main = "Variogramme", xlab = "Distance", ylab = "Semi-variance")
points(emp_vario$dist, emp_vario$vario, pch = 16, col = "red", cex = 1.2)
legend("bottomright", c("Modèle", "Empirique"), 
       col = c("blue", "red"), lty = c(1, NA), pch = c(NA, 16), 
       lwd = 2, cex = 0.8, bty = "n")

# 4. Carte des sites avec distances
plot(coord[, 1], coord[, 2], 
     pch = 16, col = "darkgreen", cex = 0.6,
     main = paste("Site locations (n =", ncol(D_Mat), ")"),
     xlab = "Longitude", ylab = "Latitude")
# Ajouter quelques liens pour montrer les distances
if(ncol(D_Mat) <= 100) {
  set.seed(123)
  sample_pairs <- sample(1:ncol(D_Mat), min(20, ncol(D_Mat)))
  for(i in 1:(length(sample_pairs)-1)) {
    segments(coord[sample_pairs[i], 1], coord[sample_pairs[i], 2],
             coord[sample_pairs[i+1], 1], coord[sample_pairs[i+1], 2],
             col = rgb(0, 0, 1, 0.1))
  }
}

# 5. Distribution des distances
hist(as.vector(dist(coord)), breaks = 50, col = "lightblue", border = "white",
     main = "Distribution of inter-site distances",
     xlab = "Distance", ylab = "Frequency")
if(exists("h_min")) {
  abline(v = c(h_min, h_max), col = "darkgreen", lwd = 2, lty = 2)
  text(h_min, par("usr")[4] * 0.1, 
       sprintf("h_min=%.2f", h_min), pos = 4, col = "darkgreen")
  text(h_max, par("usr")[4] * 0.9, 
       sprintf("h_max=%.2f", h_max), pos = 2, col = "darkgreen")
}

# 6. Paramètres du modèle
param_names <- names(extt$param)
param_values <- extt$param
barplot(param_values, names.arg = param_names, 
        col = "steelblue", las = 2,
        main = "Estimated model parameters",
        ylab = "Value", cex.names = 0.8)

par(mfrow = c(1, 1))

dev.off()

# Sauvegarder aussi en PNG
png("diagnostic_maxstab_complet.png", width = 16, height = 10, 
    units = "in", res = 300)

par(mfrow = c(2, 3), mar = c(4.5, 4.5, 3, 2))

# Répéter les mêmes graphiques...
# (même code que ci-dessus)

dev.off()

cat("Graphiques sauvegardés: diagnostic_maxstab_complet.pdf et .png\n\n")

## 3.9. TABLEAU RÉCAPITULATIF---------

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                  TABLEAU RÉCAPITULATIF                        ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

recap <- data.frame(
  Critère = c("CLIC", "CLICc", "BIC", 
              "RMSE (Ext. Coef.)", "MAE (Ext. Coef.)", "Corr (Ext. Coef.)",
              "biais", "Convergence"),
  Valeur = c(
    round(clic_results$CLIC, 2),
    round(clic_results$CLICc, 2),
    round(clic_results$BIC, 2),
    round(rmse_ext, 4),
    round(mae_ext, 4),
    round(cor_ext, 4),
    round(biais_ext, 4),
    extt_t$convergence
  )
)

print(recap, row.names = FALSE)

cat("\n")

# Tableau des distances pour θ(h)
if(exists("h_1.3")) {
  cat("Distances pour les seuils de θ(h):\n")
  tableau_h <- data.frame(
    `θ(h)` = c(1.3, 1.5, 1.7),
    Distance_h = c(h_1.3, h_1.5, h_1.7)
  )
  print(tableau_h, row.names = FALSE)
}

cat("\n")

## 3.10. EXPORT DES RÉSULTATS

# Sauvegarder les résultats numériques
resultats_export <- list(
  info_modele = list(
    type = "T-Power Exponential",
    n_sites = ncol(D_Mat),
    n_annees = nrow(D_Mat),
    convergence = extt_t$convergence
  ),
  parametres = extt_t$param,
  erreurs_standard = if(!is.null(extt_t$std.err)) extt_t$std.err else NULL,
  criteres = clic_results,
  metriques_validation = list(
    RMSE_ext = rmse_ext,
    MAE_ext = mae_ext,
    biais_ext = biais_ext,
    Correlation_ext = cor_ext
  ),
  distances_theta = if(exists("h_1.3")) {
    list(
      h_1.3 = h_1.3,
      h_1.5 = h_1.5,
      h_1.7 = h_1.7,
      h_min = h_min,
      h_max = h_max
    )
  } else NULL,
  coef_extremal_theorique = df_theorique,
  coef_extremal_empirique = data.frame(
    distance = emp_ext[,1],
    theta = emp_ext[,3]
  )
)

saveRDS(resultats_export, "resultats_diagnostic_maxstab_extt.rds")

# Export CSV
write.csv(df_theorique, "coefficient_extremal_theorique.csv", row.names = FALSE)
write.csv(data.frame(distance = emp_ext[,1], theta = emp_ext[,3]),
          "coefficient_extremal_empirique.csv", row.names = FALSE)

if(exists("h_1.3")) {
  write.csv(tableau_h, "distances_theta_seuils.csv", row.names = FALSE)
}

cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║            DIAGNOSTIC COMPLET TERMINÉ                         ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("Fichiers générés:\n")
cat("  - diagnostic_maxstab_complet.pdf\n")
cat("  - diagnostic_maxstab_complet.png\n")
cat("  - resultats_diagnostic_maxstab.rds\n")
cat("  - coefficient_extremal_theorique.csv\n")
cat("  - coefficient_extremal_empirique.csv\n")
if(exists("h_1.3")) {
  cat("  - distances_theta_seuils.csv\n")
}

cat("\n✓ Analyse terminée avec succès!\n")








#4. Cas train test---------
# ==
## 4.1 : AJUSTEMENT DU MODÈLE SUR LES 1255 SITES----
# ==
# Extraire les coordonnées pour les deux sous-ensembles
coord_sample <- coord[stations_sub, ]
coord_reste <- coord[stations_reste, ]

# Matrice temporelle pour l'échantillon
temp.cov_sample <- temp.cov[, stations_sub]

# Ajuster le modèle max-stable sur les 1255 sites
start_sample <- c(list(range = 2, smooth = 0.5, DoF = 7.5), 
                  Mod3$fitted.values[1:9])  # Adapter selon vos paramètres

extt_train_fit <- fitmaxstab(D_Mat_sample, coord_sample, "tpowexp", 
                             loc.form_3, scale.form_3, shape.form_3, 
                             temp.cov = temp.cov_sample, 
                             temp.form.loc = temp.form.loc_3, 
                             temp.form.scale = temp.form.scale_3,
                             nugget = 0, method = "nlm",
                             start = start_sample)

# ==
## 4.2 : INTERPOLATION SUR LES 313 SITES RESTANTS-------
# ==
# Préparer les nouvelles données pour la prédiction
# Il faut créer un data.frame avec les coordonnées et les covariables temporelles

# Créer un data.frame pour chaque combinaison site-année
newdata_reste <- data.frame(
  lat = rep(coord_reste[, 1]),# each = 24),
  lon = rep(coord_reste[, 2])#, each = 24),
  #years = rep(1:24, times = length(stations_reste))
)

# Prédiction des paramètres GEV et niveau de retour à 50 ans
pred_results <- predict(extt_train_fit, 
                        newdata = newdata_reste, 
                        ret.per = 50,  # Niveau de retour à 50 ans
                        std.err = TRUE)
#View(pred_results)
pred_results = as.data.frame(pred_results)

## Modèle gev simple sur data reste
m = c(); sc = c();  sh = c()
for(i in seq(1, ncol(D_Mat_reste), 1)){
  param = gevmle(D_Mat_reste[,i])
  m[i] = as.numeric(param[1])
  sc[i] = as.numeric(param[2])
  sh[i] = as.numeric(param[3])
}

df_reste = data.frame(lat = pred_results$lat,
                      lon = pred_results$lon,
                      mu = m,
                      sigma = sc,
                      xi = sh)



# Extraire les valeurs prédites
# La structure exacte dépend du package, typiquement:
# pred_results$fit contient les prédictions
# pred_results$se.fit contient les erreurs standard

# Organiser les prédictions en matrice (24 années x 313 sites)
loc_pred <- matrix(pred_results$loc, nrow = 24, ncol = length(stations_reste))
scale_pred <- matrix(pred_results$scale, nrow = 24, ncol = length(stations_reste))
shape_pred <- matrix(pred_results$shape, nrow = 24, ncol = length(stations_reste))

# Si disponible, extraire les niveaux de retour prédits
if (!is.null(pred_results$Q50)) {
  niveaux_retour_pred <- matrix(pred_results$Q50, 
                                nrow = 24, 
                                ncol = length(stations_reste))
  # Prendre la dernière année (2023)
  niveaux_retour_pred_2023 <- niveaux_retour_pred[24, ]
}

# ==
## 4.3 : CALCUL DES VALEURS PRÉDITES (MÉDIANES GEV) ----
# ==
# Fonction pour calculer la médiane d'une GEV
median_gev <- function(loc, scale, shape) {
  if (abs(shape) < 1e-6) {
    return(loc - scale * log(-log(1-1/2)))
  } else {
    return(loc + (scale/shape) * ((-log(1-1/2))^(-shape) - 1))
  }
}

# Calculer les prédictions (médiane) pour chaque site et année
D_Mat_pred <- matrix(NA, nrow = 24, ncol = length(stations_reste))

for (i in 1:ncol(D_Mat_pred)) {
  for (t in 1:24) {
    D_Mat_pred[t, i] <- median_gev(loc_pred[t, i], 
                                   scale_pred[t, i], 
                                   shape_pred[t, i])
  }
}

# ==
## 4.4 : CALCUL DES MÉTRIQUES DE PERFORMANCE-----
# ==
library(dplyr)

# Observations réelles vs prédictions
obs <- as.vector(D_Mat_reste)
pred <- as.vector(D_Mat_pred)

# Retirer les NA
valid_idx <- !is.na(obs) & !is.na(pred) & is.finite(obs) & is.finite(pred)
obs_clean <- obs[valid_idx]
pred_clean <- pred[valid_idx]

# Métriques globales
metriques_globales <- data.frame(
  RMSE = sqrt(mean((obs_clean - pred_clean)^2)),
  MAE = mean(abs(obs_clean - pred_clean)),
  Biais = mean(pred_clean - obs_clean),
  Biais_relatif_pct = mean((pred_clean - obs_clean) / obs_clean) * 100,
  Correlation = cor(obs_clean, pred_clean),
  R2 = cor(obs_clean, pred_clean)^2,
  NSE = 1 - sum((obs_clean - pred_clean)^2) / sum((obs_clean - mean(obs_clean))^2)
)

cat("\n=== MÉTRIQUES DE PERFORMANCE GLOBALES ===\n")
print(round(metriques_globales, 4))

# Métriques par station
df_validation <- data.frame(
  station = rep(stations_reste, each = 24),
  annee = rep(1:24, length(stations_reste)),
  obs = as.vector(D_Mat_reste),
  pred = as.vector(D_Mat_pred)
) %>%
  filter(!is.na(obs) & !is.na(pred) & is.finite(obs) & is.finite(pred)) %>%
  mutate(
    erreur = pred - obs,
    erreur_abs = abs(erreur),
    erreur_carre = erreur^2,
    erreur_relative_pct = (erreur / obs) * 100
  )

metriques_par_station <- df_validation %>%
  group_by(station) %>%
  summarise(
    n_obs = n(),
    RMSE = sqrt(mean(erreur_carre)),
    MAE = mean(erreur_abs),
    Biais = mean(erreur),
    Biais_relatif_pct = mean(erreur_relative_pct),
    Correlation = cor(obs, pred),
    R2 = Correlation^2,
    .groups = 'drop'
  )

cat("\n=== STATISTIQUES DES MÉTRIQUES PAR STATION ===\n")
print(summary(metriques_par_station[, -1]))

# Métriques par année
metriques_par_annee <- df_validation %>%
  group_by(annee) %>%
  summarise(
    n_obs = n(),
    RMSE = sqrt(mean(erreur_carre)),
    MAE = mean(erreur_abs),
    Biais = mean(erreur),
    Correlation = cor(obs, pred),
    .groups = 'drop'
  )

cat("\n=== MÉTRIQUES PAR ANNÉE ===\n")
print(metriques_par_annee)

# ==
##  4.5 : NIVEAUX DE RETOUR À 50 ANS ------
# ==
library(evd)

# Fonction pour calculer le niveau de retour
calcul_niveau_retour <- function(loc, scale, shape, periode = 50) {
  if (abs(shape) < 1e-6) {
    return(loc - scale * log(-log(1 - 1/periode)))
  } else {
    return(loc + (scale/shape) * ((-log(1 - 1/periode))^(-shape) - 1))
  }
}

# Si les niveaux de retour n'ont pas été calculés par predict()
if (!exists("niveaux_retour_pred_2023")) {
  niveaux_retour_pred_2023 <- sapply(1:ncol(D_Mat_pred), function(i) {
    calcul_niveau_retour(loc_pred[24, i], scale_pred[24, i], shape_pred[24, i], 50)
  })
}

# Calculer les niveaux de retour observés
niveaux_retour_obs <- sapply(1:length(stations_reste), function(i) {
  data_station <- D_Mat_reste[, i]
  
  tryCatch({
    fit <- fgev(data_station[!is.na(data_station)])
    calcul_niveau_retour(fit$estimate[1], fit$estimate[2], fit$estimate[3], 50)
  }, error = function(e) {
    return(NA)
  })
})

# Comparaison des niveaux de retour
comparaison_retour <- data.frame(
  station = stations_reste,
  lat = coord_reste[, 1],
  lon = coord_reste[, 2],
  niveau_retour_obs = niveaux_retour_obs,
  niveau_retour_pred = niveaux_retour_pred_2023
) %>%
  filter(!is.na(niveau_retour_obs) & !is.na(niveau_retour_pred)) %>%
  mutate(
    erreur = niveau_retour_pred - niveau_retour_obs,
    erreur_abs = abs(erreur),
    erreur_relative_pct = (erreur / niveau_retour_obs) * 100
  )

cat("\n=== STATISTIQUES SUR LES NIVEAUX DE RETOUR À 50 ANS ===\n")
print(summary(comparaison_retour[, c("niveau_retour_obs", "niveau_retour_pred", 
                                     "erreur", "erreur_relative_pct")]))

metriques_retour <- data.frame(
  RMSE = sqrt(mean(comparaison_retour$erreur^2)),
  MAE = mean(comparaison_retour$erreur_abs),
  Biais = mean(comparaison_retour$erreur),
  Biais_relatif_pct = mean(comparaison_retour$erreur_relative_pct),
  Correlation = cor(comparaison_retour$niveau_retour_obs, 
                    comparaison_retour$niveau_retour_pred)
)

cat("\n=== MÉTRIQUES POUR LES NIVEAUX DE RETOUR ===\n")
print(round(metriques_retour, 4))

# ==
## 4.6 : VISUALISATIONS ------
# ==
library(ggplot2)
library(gridExtra)

# 1. Observations vs Prédictions (toutes années)
p1 <- ggplot(df_validation, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen", alpha = 0.2) +
  labs(title = "Observations vs Prédictions (313 sites, 24 ans)",
       subtitle = paste0("R² = ", round(metriques_globales$R2, 3), 
                         " | RMSE = ", round(metriques_globales$RMSE, 2)),
       x = "Observations (mm)", 
       y = "Prédictions (mm)") +
  theme_minimal(base_size = 12)

# 2. Distribution des erreurs
p2 <- ggplot(df_validation, aes(x = erreur)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = median(df_validation$erreur), 
             color = "darkgreen", linetype = "dotted", size = 1) +
  labs(title = "Distribution des erreurs de prédiction",
       subtitle = paste0("Biais = ", round(metriques_globales$Biais, 2), " mm"),
       x = "Erreur (Pred - Obs) en mm", 
       y = "Fréquence") +
  theme_minimal(base_size = 12)

# 3. RMSE par station (carte ou barplot)
p3 <- ggplot(metriques_par_station %>% 
               arrange(desc(RMSE)) %>% 
               head(30), 
             aes(x = reorder(factor(station), RMSE), y = RMSE)) +
  geom_col(fill = "coral", alpha = 0.8) +
  coord_flip() +
  labs(title = "Top 30 des stations avec le plus grand RMSE",
       x = "Station", 
       y = "RMSE (mm)") +
  theme_minimal(base_size = 10)

# 4. Évolution temporelle du biais
p4 <- ggplot(metriques_par_annee, aes(x = annee, y = Biais)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Évolution du biais au cours du temps",
       x = "Année", 
       y = "Biais (mm)") +
  theme_minimal(base_size = 12)

# 5. Niveaux de retour: Obs vs Pred
p5 <- ggplot(comparaison_retour, aes(x = niveau_retour_obs, y = niveau_retour_pred)) +
  geom_point(alpha = 0.6, color = "darkgreen", size = 2.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.2) +
  labs(title = "Niveaux de retour à 50 ans: Obs vs Pred",
       subtitle = paste0("R² = ", round(metriques_retour$Correlation^2, 3), 
                         " | RMSE = ", round(metriques_retour$RMSE, 2)),
       x = "Niveau de retour observé (mm)", 
       y = "Niveau de retour prédit (mm)") +
  theme_minimal(base_size = 12)

# 6. Erreurs relatives sur niveaux de retour
p6 <- ggplot(comparaison_retour, aes(x = erreur_relative_pct)) +
  geom_histogram(bins = 40, fill = "darkgreen", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = median(comparaison_retour$erreur_relative_pct), 
             color = "blue", linetype = "dotted", size = 1) +
  labs(title = "Erreurs relatives sur les niveaux de retour",
       x = "Erreur relative (%)", 
       y = "Fréquence") +
  theme_minimal(base_size = 12)

# Afficher les graphiques
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
print(p6)

# Graphique combiné
grid.arrange(p1, p2, p4, p5, ncol = 2)

# ==
## 4.7 : SAUVEGARDE DES RÉSULTATS----
# ==
# Sauvegarder les résultats
write.csv(metriques_globales, "resultats/metriques_globales.csv", row.names = FALSE)
write.csv(metriques_par_station, "resultats/metriques_par_station.csv", row.names = FALSE)
write.csv(metriques_par_annee, "resultats/metriques_par_annee.csv", row.names = FALSE)
write.csv(comparaison_retour, "resultats/comparaison_niveaux_retour_50ans.csv", row.names = FALSE)
write.csv(df_validation, "resultats/donnees_validation_complete.csv", row.names = FALSE)

# Sauvegarder le modèle
saveRDS(extt_t_sample, "resultats/modele_maxstab_1255sites.rds")

cat("\n=== ANALYSE TERMINÉE ===\n")
cat("Tous les résultats ont été sauvegardés dans le dossier 'resultats/'\n")

# 5. Niveau de retour cas dépendant--------
## a. Toutes les années --------
Niveau_Retour= function(la, lo, t, Annee){
  time=0 ; z_T = 0
  if (Annee <=2023 & Annee >=2000){
    time =(Annee-2000+1)
    
    #Meilleur modèle spatial
    mu = extt_loc_intercept + extt_loc_lat*la + extt_loc_lon*lo 
    sigma = extt_scale_intercept + extt_scale_lat*la + extt_scale_lon*lo 
    xi = extt_shape_intercept 
    
    
    #Meilleur modèle spatio temporel sans dépendance sinusoïdale et xi dpt uniqmt du tps
    #mu = Mod_1_loc_intercept + Mod_1_loc_lat*la + Mod_1_loc_lon*lo + Mod_1_temp_loc*time 
    #sigma = Mod_1_scale_intercept + Mod_1_scale_lat*la + Mod_1_scale_lon*lo + Mod_1_temp_scale*time
    #xi = Mod_1_shape_intercept
    
    #Meilleur modèle spatio temporel sans dépendance sinusoïdale
    #mu = Mod1_loc_intercept + Mod1_loc_lat*la + Mod1_loc_lon*lo + Mod1_temp_loc*time 
    #sigma = Mod1_scale_intercept + Mod1_scale_lat*la + Mod1_scale_lon*lo + Mod1_temp_scale*time
    #xi = Mod1_shape_intercept 
    
    
    #Meilleur modèle spatio temporel avec periode 24
    #mu = M6_loc_intercept + M6_loc_lat*la + M6_loc_lon*lo + M6_temp_loc*time 
    #sigma = M6_scale_intercept + M6_scale_lat*la + M6_scale_lon*lo + M6_temp_scale*time
    #xi = M6_shape_intercept +  M6_temp_shape_sin*time
    
    #Meilleur modèle spatio temporel avec periode 12
    #mu = M8_loc_intercept + M8_loc_lat*la + M8_loc_lon*lo + M8_temp_loc*time 
    #sigma = M8_scale_intercept + M8_scale_lat*la + M8_scale_lon*lo + M8_temp_scale*time
    #xi = M8_shape_intercept +  M8_temp_shape_cos*time
    
    #Meilleur modèle spatio temporel avec periode 5
    #mu = M1_loc_intercept + M1_loc_lat*la + M1_loc_lon*lo + M1_temp_loc*time 
    #sigma = M1_scale_intercept + M1_scale_lat*la + M1_scale_lon*lo + M1_temp_scale*time
    #xi = M1_shape_intercept +  M1_temp_shape*time
    
    
    #Meilleur modèle spatio temporel avec normalisation de lat et lon et periode 24
    #mu = M_sc_loc_intercept + M_sc_loc_lat*la + M_sc_loc_lon*lo + M_sc_temp_loc*time 
    #sigma = M_sc_scale_intercept + M_sc_scale_lat*la + M_sc_scale_lon*lo + M_sc_temp_scale*time
    #xi = M_sc_shape_intercept +  M_sc_temp_shape_sin*time
    
    if (abs(xi) <= 0.001) {
      # Cas Gumbel : xi = 0
      z_T <- mu - sigma * log(-log(1 - 1 / t))
    } else {
      # Cas général : xi != 0
      z_T <- mu + (sigma/xi) * ( (-log(1 - 1 / t))^(-xi) - 1)
    }
  }
  
  list(MU=mu, Sig= sigma, Rtn =z_T)
}

AN = 2000:2023

COORD <- matrix( c(PL$LAT,PL$LON), ncol=2) 
COORD_mult= COORD[rep(1:nrow(COORD), 24), ]

Result = function(t){
  Mu_m = Sig_m= Return = matrix(0, 1568*24 ,4)
  for(k in 1:24){
    for(i in 1:1568){
      Res = Niveau_Retour(la=COORD[i,1], lo=COORD[i,2], t=t, Annee=AN[k])
      
      Mu_m[1568*(k-1)+i,] = c(COORD[i,], Res$MU, AN[k])
      Sig_m[1568*(k-1)+i,] = c(COORD[i,], Res$Sig, AN[k])
      Return[1568*(k-1)+i,]  = c(COORD[i,], Res$Rtn, AN[k])
    }
  }
  list(Mu_mat =Mu_m , Sig_mat=Sig_m, Return.level=Return)
}

## b. Dans T ans-------
compute_return_levels <- function(la, lo, T = c(10,50,100), tol_xi = 1e-8){
  
  #Meilleur modèle spatial
  #mu = extt_loc_intercept + extt_loc_lat*la + extt_loc_lon*lo 
  #sigma = extt_scale_intercept + extt_scale_lat*la + extt_scale_lon*lo 
  #xi = extt_shape_intercept 
  
  mu = extt1_loc_intercept + extt1_loc_lat*la + extt1_loc_lon*lo 
  sigma = extt1_scale_intercept + extt1_scale_lat*la + extt1_scale_lon*lo 
  xi = extt1_shape_intercept 

  if(any(sigma <= 0)){
    warning("Certaines valeurs de sigma <= 0 : vérifier la parametrisation ou contraindre sigma>0.")
  }
  
  q <- -log(1 - 1 / T)  # positive numbers, length = length(T)
  
  # output: array (n_sites x n_T)
  la_vec <- as.vector(la); n_sites <- length(la_vec)
  res <- matrix(NA_real_, nrow = n_sites, ncol = length(T),
                dimnames = list(NULL, paste0("T", T)))
  
  for(iT in seq_along(T)){
    qi <- q[iT]
    if(abs(xi) > tol_xi){
      # general case
      res[, iT] <- mu + (sigma / xi) * ( qi^(-xi) - 1 )
    } else {
      # xi ~ 0 : Gumbel limit
      res[, iT] <- mu - sigma * log(qi)
    }
    # check domain 1 + xi*(z-mu)/sigma > 0 (optional)
    # but formula already ensures it if sigma>0 and q>0
  }
  as.data.frame(res)
}

ReturnLevel = compute_return_levels(precip$Latitude, precip$Longitude, T = c(10,50,100), tol_xi = 1e-8)
ReturnLevel$Latitude = precip$Latitude
ReturnLevel$Longitude = precip$Longitude
#write.csv2(ReturnLevel, "extt_RL2.csv")
write.csv2(ReturnLevel, "extt1_RL.csv")
# 6. Carte des niveaux de retour---------
## a Cas1---------
RESULTAT10=Result(10)
RESULTAT50=Result(50)
RESULTAT100=Result(100)

Niveaux_retour  <- data.frame(
  Latitude = RESULTAT100$Mu_mat[,1],
  Longitude = RESULTAT100$Mu_mat[,2],
  Mu = RESULTAT100$Mu_mat[,3], 
  Sig = RESULTAT100$Sig_mat[,3],
  Year  =  RESULTAT100$Mu_mat[,4],
  Return_level10 = RESULTAT10$Return.level[,3],
  Return_level50 = RESULTAT50$Return.level[,3],
  Return_level100 = RESULTAT100$Return.level[,3]
)
attach(Niveaux_retour)

#write.csv2(Niveaux_retour, 'extt_RL.csv')
#View(Niveaux_retour)

library(RColorBrewer) ; library(ggplot2) ; library(latex2exp); library(maps) 

pays = c("Benin", "Burkina Faso", "Cabo Verde", "Ivory Coast", 
         "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
         "Liberia", "Mali", "Mauritania", "Niger", "Nigeria", 
         "Senegal", "Sierra Leone", "Togo")

ggplot() +
  # Couleur des frontières
  geom_point(data = Niveaux_retour, aes(x = Longitude, y = Latitude, color =  Return_level50),
             size = 8, stat = "identity") + 
  #scale_color_gradient(low = c("darkblue","lightblue"), high = c("red")) + # Dégradé de couleur
  
  scale_color_gradientn(colors = c("lightblue","darkblue", "green", "yellow", "orange", "red")) + # Dégradé de couleur
  #scale_fill_viridis_d(option = "A")+
  
  scale_x_continuous(name = "Longitude", 
                     breaks = seq(-17, 11, by = 6),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(name = "Latitude", 
                     breaks = seq(4, 20, by = 3),
                     labels = function(y) paste0(y, "°N")) +
  
  labs(title = "50-year return level", 
       #  x = "Longitude",
       #    y = "Latitude",
       color = "mm/jr")+
  theme_minimal() +
  #facet_wrap(~ Year) +
  
  borders("world", regions = pays, colour = "black",size=.8) +
  
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

## b Cas2-----------

ggplot() +
  # Couleur des frontières
  geom_point(data = ReturnLevel, aes(x = Longitude, y = Latitude, color =  T100), size = 8, stat = "identity") + 
  #scale_color_gradient(low = c("darkblue","lightblue"), high = c("red")) + # Dégradé de couleur
  
  scale_color_gradientn(colors = c("lightblue","darkblue", "green", "yellow", "orange", "red")) + # Dégradé de couleur
  #scale_fill_viridis_d(option = "A")+
  
  scale_x_continuous(name = "Longitude", 
                     breaks = seq(-17, 11, by = 6),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(name = "Latitude", 
                     breaks = seq(4, 20, by = 3),
                     labels = function(y) paste0(y, "°N")) +
  
  labs(title = "100-year return level", 
       #  x = "Longitude",
       #    y = "Latitude",
       color = "Level(mm/jr)")+
  theme_minimal() +
  
  borders("world", regions = c("Benin", "Burkina Faso", "Cabo Verde", "Ivory Coast", 
                               "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
                               "Liberia", "Mali", "Mauritania", "Niger", "Nigeria", 
                               "Senegal", "Sierra Leone", "Togo"),
          colour = "black",size=.8) +
  
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


save.image(file = "maxstable")