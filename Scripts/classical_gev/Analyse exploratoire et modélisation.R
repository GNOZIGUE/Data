load("analyse_exploratoire.rds")

library(SpatialExtremes); library(RColorBrewer); library(ggplot2); library(latex2exp)
library(readxl); library(tidyr); library(dplyr); library(patchwork)
precipitation <- read_excel("../Matrice_des_parametres.xlsx", 
                            sheet = "Pr\u00E9cipitation")
#View(precipitation)
precipTr <- read_excel("../Matrice_des_parametres.xlsx", 
                       sheet = "Transpos\u00E9Precipitation")
#View(precipTr)

PL <- read.delim("../AO.txt")
#View(PL)
attach(PL)

coord <- matrix(c(LAT,LON), ncol=2) ;  colnames(coord) <- c("lat", "lon")

D_Mat = matrix(c(A2000, A2001, A2002, A2003, A2004, A2005, A2006, A2007, A2008, A2009,
                 A2010, A2011, A2012, A2013, A2014 ,A2015 ,A2016, A2017, A2018, A2019, A2020, A2021, A2022 , A2023),
               ncol=1568, byrow=TRUE )
#View(D_Mat)

# Série chronologique -------

data_long <- precipitation %>%
  pivot_longer(cols = starts_with("20"), 
               names_to = "Year", 
               values_to = "Precipitation")
#View(data_long)
attach(data_long)



ggplot(data_long, aes(x = Year, y = Precipitation, colour = Precipitation)) +
  geom_bar(stat = "identity") +
  #scale_y_continuous(breaks = seq(0, 1000, 200)) +
  scale_x_discrete(breaks = seq(2000, 2023, 5)) +
  labs(title = "Evolution temporelle de la pluviométrie", 
       x = "Années",
       y = "Précipitation")+
  theme_minimal() + 
  facet_wrap(~ LAT) +
  theme(legend.title = element_text(),
        axis.text.x = element_text(size = 8, face = "bold", angle = 90),
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 10, face = "bold"))

#A- Analyse exploratoire spatiale -------

# 1) Distribution spatiale des paramètres gev--------

## Modèle gev simple 
mu = c(); sig = c();  psi = c()
for(i in seq(1, ncol(D_Mat), 1)){
  param = gevmle(D_Mat[,i])
  mu[i] = as.numeric(param[1])
  sig[i] = as.numeric(param[2])
  psi[i] = as.numeric(param[3])
}

precip_st = data.frame(Latitude = precipitation$LAT,
                       Longitude = precipitation$LON,
                       Location = mu, 
                       Scale = sig,
                       shape = psi)

#View(precip_st)
attach(precip_st)


fig_loc = ggplot(precip_st) +
  aes(x = Longitude, y = Latitude, colour = Location) +
  geom_point(size = 5, shape = "square") +
  scale_color_distiller(palette = "RdYlGn", direction = -1) +
  scale_x_continuous(breaks = seq(-17, 11, by = 7),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(breaks = seq(4, 20, by = 2),
                     labels = function(y) paste0(y, "°N")) +
  labs(x = "Longitude", y = "Latitude", 
       title = "Location parameter") +
  theme_minimal() +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.position="right", 
        legend.text = element_text(size = 8, face = "bold"),
        #legend.position.inside =c(15,18), 
        plot.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"), 
        axis.title.y = element_text(size = 10, face = "bold"))+
  annotation_borders("world", regions = c("Benin", "Burkina Faso", "Ivory Coast", 
                               "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
                               "Liberia", "Mali", "Mauritania", "Niger", 
                               "Nigeria", "Senegal", "Sierra Leone", "Togo"), 
          colour = "black",size=.8) +
  coord_sf(ylim = c(3, 18), xlim = c(-17, 10.5)) 

###


fig_loc2 = ggplot() +
  # Couleur des frontières
  geom_point(data = precip_st, aes(x = Longitude, y = Latitude, color =  Location),
             size =15,stat = "identity", shape = "square") + # Points colorés par température
  #scale_color_gradient(low = c("darkblue","lightblue"), high = c("red")) + # Dégradé de couleur
  
  scale_color_gradient(low = c("darkblue","lightblue", "green"), high = c("yellow", "orange", "red")) + # Dégradé de couleur
  #scale_fill_viridis_d(option = "A")+
  
  scale_x_continuous(name = "Longitude", 
                     breaks = seq(-17, 11, by = 3),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(name = "Latitude", 
                     breaks = seq(4, 20, by = 3),
                     labels = function(y) paste0(y, "°N")) +
  
  labs(title = "Location parameter", 
       #  x = "Longitude",
       #    y = "Latitude",
       color = "Value")+
  theme_minimal() +
  
  borders("world", regions = c("Benin", "Burkina Faso", "Ivory Coast", 
                               "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
                               "Liberia", "Mali", "Mauritania", "Niger", 
                               "Nigeria", "Senegal", "Sierra Leone", "Togo"), 
          colour = "black",size=.8) +
  
  theme(legend.title = element_text(hjust = 0.5, size = 17,
                                    face = "bold"), 
        legend.position="right",
        #legend.position.inside =c(15,18), 
        plot.title = element_text(hjust = 0.5, size = 18),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 17, face = "bold"), 
        axis.title.y = element_text(size = 17, face = "bold"))+
  coord_fixed( xlim = c(-16.75, 10.1), ylim = c(4.25, 17.75))  


########################-

fig_sig = ggplot(precip_st) +
  aes(x = Longitude, y = Latitude, colour = Scale) +
  geom_point(size = 5, shape = "square") +
  scale_color_distiller(palette = "RdYlGn", direction = -1) +
  scale_x_continuous(breaks = seq(-17, 11, by = 7),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(breaks = seq(4, 20, by = 2),
                     labels = function(y) paste0(y, "°N")) +
  labs(x = "Longitude", y = "Latitude", 
       title = "Scale parameter") +
  theme_minimal() +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.position="right", 
        legend.text = element_text(size = 8, face = "bold"),
        #legend.position.inside =c(15,18), 
        plot.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"), 
        axis.title.y = element_text(size = 10, face = "bold"))+
  borders("world", regions = c("Benin", "Burkina Faso", "Ivory Coast", 
                               "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
                               "Liberia", "Mali", "Mauritania", "Niger", 
                               "Nigeria", "Senegal", "Sierra Leone", "Togo"), 
          colour = "black",size=.8) +
  coord_sf(ylim = c(3, 18), xlim = c(-17, 10.5)) 

fig_psi = ggplot(precip_st) +
  aes(x = Longitude, y = Latitude, colour = shape) +
  geom_point(size = 5, shape = "square") +
  scale_color_distiller(palette = "RdYlGn", direction = -1) +
  scale_x_continuous(breaks = seq(-17, 11, by = 7),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(breaks = seq(4, 20, by = 2),
                     labels = function(y) paste0(y, "°N")) +
  labs(x = "Longitude", y = "Latitude", 
       title = "Shape parameter") +
  theme_minimal() +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.position="right", 
        legend.text = element_text(size = 8, face = "bold"),
        #legend.position.inside =c(15,18), 
        plot.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 8, face = "bold"),
        axis.title.x = element_text(size = 10, face = "bold"), 
        axis.title.y = element_text(size = 10, face = "bold"))+
  borders("world", regions = c("Benin", "Burkina Faso", "Cote-d'Ivoire", 
                               "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
                               "Liberia", "Mali", "Mauritania", "Niger", 
                               "Nigeria", "Senegal", "Sierra Leone", "Togo"), 
          colour = "black",size=.8) +
  coord_sf(ylim = c(3, 18), xlim = c(-17, 10.5)) 


#fig_loc + fig_sig + fig_psi 
fig_loc / fig_sig / fig_psi

#2) Paramètres estimés en fonction de la longitude ---------

loc_fig = ggplot(precip_st) +
  aes(x = Longitude, y = Location) +
  geom_jitter() +
  geom_smooth(method = "glm", formula = y ~ x) +
  scale_x_continuous(breaks = seq(-17, 11, by = 2),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  labs(x = "Longitude", y = TeX("$mu$"), 
       title = "Location") +
  theme_minimal()

sig_fig = ggplot(precip_st) +
  aes(x = Longitude, y = Scale) +
  geom_jitter() +
  geom_smooth(method = "glm", formula = y ~ x) +
  scale_x_continuous(breaks = seq(-17, 11, by = 2),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  labs(x = "Longitude", y = TeX("$sigma$"), 
       title = "Scale") +
  theme_minimal()

psi_fig = ggplot(precip_st) +
  aes(x = Longitude, y = shape) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  scale_x_continuous(breaks = seq(-17, 11, by = 2),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  labs(x = "Longitude", y = TeX("$xi$"), 
       title = "Shape") +
  theme_minimal()

loc_fig / sig_fig / psi_fig

#3) Paramètres estimés en fonction de la latitude ---------
fig_loc = ggplot(precip_st) +
  aes(x = Latitude, y = Location) +
  # geom_point(size = 10, shape = "square") +
  geom_jitter() +
  geom_smooth(method = "glm", formula = y ~ x) +
  #scale_color_distiller(palette = "RdYlGn", direction = -1) +
  #scale_x_continuous(breaks = seq(-17, 11, by = 7),
  #           labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_x_continuous(breaks = seq(4, 20, by = 2),
                     labels = function(y) paste0(y, "°N")) +
  labs(x = "Latitude", y = TeX("$mu$"), 
       title = "Location") +
  theme_minimal() 

fig_sig = ggplot(precip_st) +
  aes(x = Latitude, y = Scale) +
  geom_jitter() +
  geom_smooth(method = "glm", formula = y ~ x) +
  #geom_point(size = 10, shape = "square") +
  #scale_color_distiller(palette = "RdYlGn", direction = -1) +
  #scale_x_continuous(breaks = seq(-17, 11, by = 7),
  #                  labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_x_continuous(breaks = seq(4, 20, by = 2),
                     labels = function(y) paste0(y, "°N")) +
  labs(x = "Latitude", y = TeX("$sigma$"), 
       title = "Scale") +
  theme_minimal() 

fig_psi = ggplot(precip_st) +
  aes(x = Latitude, y = shape) +
  geom_jitter() +
  geom_smooth(method = "glm", formula = y ~ x) +
  #scale_x_continuous(breaks = seq(-17, 11, by = 7),labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_x_continuous(breaks = seq(4, 20, by = 2), 
                     labels = function(y) paste0(y, "°N")) +
  labs(x = "Latitude", y = TeX("$xi$"), 
       title = "Shape") +
  theme_minimal()

fig_loc / fig_sig / fig_psi


#4) Dépendance des paramètres GEV: ----------
###a) En fonction de la latitude -------
# On estime ces paramètres pour chaque valeur de latitude en chaque année 

Bande_Lat = seq(4.25, 17.75, 0.5)

precip_lat = split(precipitation[,-c(1,2,3)], (seq(nrow(precipitation)) - 1) %/% 56)

# Initialiser un vecteur pour stocker les résultats

mu_lat = numeric(length(precip_lat))
sig_lat = numeric(length(precip_lat))
psi_lat = numeric(length(precip_lat))
# Ces vecteurs contiennent des zéros qu'il faut supprimer après

# Boucle sur chaque bloc
for (i in seq_along(precip_lat)) {
  bloc <- precip_lat[[i]]
  
  # Boucle sur chaque colonne du bloc
  for (colonne in colnames(bloc)) {
    # Extraire la colonne
    valeurs <- bloc[[colonne]]
    
    #gev suivant mle sur la precip en fct de la lat
    fit = gevmle(valeurs)
    mu = as.numeric(fit[1])
    sig = as.numeric(fit[2])
    psi = as.numeric(fit[3])
    
    # Stocker les résultats dans une liste
    mu_lat <- append(mu_lat, mu)
    sig_lat <- append(sig_lat, sig)
    psi_lat <- append(psi_lat, psi)
  }
  
}

# Suppression des  valeurs initiales 0
MU_Lat = mu_lat[-c(1:28)] 
SIG_Lat = sig_lat[-c(1:28)]
PSI_Lat = psi_lat[-c(1:28)]


###### Pour la visualisation, voir Trace_mu_Sig.R


# Ordre de Mu en fct des 28 latitudes de chaque année
#MU_Lat = Mu_lat[c(seq(24,672,24), seq(23,671,24), seq(22,670,24), seq(21,669,24),
#seq(20,668,24), seq(19,667,24), seq(18,666,24), seq(17,665,24),
#seq(16,664,24), seq(15,663,24), seq(14,662,24), seq(13,661,24),
#seq(12,660,24), seq(11,659,24), seq(10,658,24), seq(9,657,24), 
#seq(8,656,24), seq(7,655,24), seq(6,654,24), seq(5,653,24),
#seq(4,652,24), seq(3,651,24), seq(2,650,24), seq(1,649,24))] 

#SIG_Lat = Sig_lat[c(seq(24,672,24), seq(23,671,24), seq(22,670,24), seq(21,669,24),
#seq(20,668,24), seq(19,667,24), seq(18,666,24), seq(17,665,24),
#seq(16,664,24), seq(15,663,24), seq(14,662,24), seq(13,661,24),
#seq(12,660,24), seq(11,659,24), seq(10,658,24), seq(9,657,24), 
#seq(8,656,24), seq(7,655,24), seq(6,654,24), seq(5,653,24),
#seq(4,652,24), seq(3,651,24), seq(2,650,24), seq(1,649,24))]

### b) En fonction de la longitude -------
#On estime ces paramètres pour chaque valeur de longitude en chaque année 
precip_lon = split(precipitation[,-c(1,2,3)], (seq(nrow(precipitation)) - 1) %/% 28)

Bande_Lon = seq(-16.75, 10.75, 0.5)

# Initialiser un vecteur pour stocker les résultats
mu_lon = numeric(length(precip_lon))
sig_lon = numeric(length(precip_lon))
psi_lon = numeric(length(precip_lon))
# Boucle sur chaque bloc
for (i in seq_along(precip_lon)) {
  bloc <- precip_lon[[i]]
  
  # Boucle sur chaque colonne du bloc
  for (colonne in colnames(bloc)) {
    # Extraire la colonne
    valeurs <- bloc[[colonne]]
    
    #gev suivant mle sur la precip en fct de la lon
    fit = gevmle(valeurs)
    mu = as.numeric(fit[1])
    sig = as.numeric(fit[2])
    psi = as.numeric(fit[3])
    
    # Stocker les résultats dans une liste
    mu_lon <- append(mu_lon, mu)
    sig_lon <- append(sig_lon, sig)
    psi_lon <- append(psi_lon, psi)
  }
  
}

# Supression des  valeurs initiales 0
MU_Lon = mu_lon[-c(1:56)] 
SIG_Lon = sig_lon[-c(1:56)]
PSI_Lon = psi_lon[-c(1:56)]


###### Pour la visualisation, voir Trace_mu_Sig.R


# Ordre de Mu en fct des 28 latitudes de chaque année
#MU_Lon = Mu_lon[c(seq(24,672,24), seq(23,671,24), seq(22,670,24), seq(21,669,24),
#seq(20,668,24), seq(19,667,24), seq(18,666,24), seq(17,665,24),
#seq(16,664,24), seq(15,663,24), seq(14,662,24), seq(13,661,24),
#seq(12,660,24), seq(11,659,24), seq(10,658,24), seq(9,657,24), 
#seq(8,656,24), seq(7,655,24), seq(6,654,24), seq(5,653,24),
#seq(4,652,24), seq(3,651,24), seq(2,650,24), seq(1,649,24))] 

#SIG_Lon = Sig_lon[c(seq(24,672,24), seq(23,671,24), seq(22,670,24), seq(21,669,24),
#seq(20,668,24), seq(19,667,24), seq(18,666,24), seq(17,665,24),
#seq(16,664,24), seq(15,663,24), seq(14,662,24), seq(13,661,24),
#seq(12,660,24), seq(11,659,24), seq(10,658,24), seq(9,657,24), 
#seq(8,656,24), seq(7,655,24), seq(6,654,24), seq(5,653,24),
#seq(4,652,24), seq(3,651,24), seq(2,650,24), seq(1,649,24))]

### c) En fonction du temps --------
annees = colnames(precipitation[,-c(1, 2, 3)])

mu_temps = numeric(length(annees))
sig_temps = numeric(length(annees))
psi_temps = numeric(length(annees))

# Boucle sur chaque colonne du bloc
for (i in seq_along(annees)) {
  # Extraire la colonne
  an <- precipitation[[annees[i]]]
  
  #gev suivant mle sur la precip en fct du temps
  fit = gevmle(an)
  mu = as.numeric(fit[1])
  sig = as.numeric(fit[2])
  psi = as.numeric(fit[3])
  
  # Stocker les résultats dans une liste
  mu_temps <- append(mu_temps, mu)
  sig_temps <- append(sig_temps, sig)
  psi_temps <- append(psi_temps, psi)
}

MU_Temps = mu_temps[-c(1:24)] 
SIG_Temps = sig_temps[-c(1:24)]
PSI_Temps = psi_temps[-c(1:24)]

#B- Modélisation spatio-temporelle --------

##1) Modélisation spatiale (longitude, latitude) ---------
#coord <- matrix(c(LAT,LON), ncol=2)

### Juste une parenthèse ##################
mul <- ~ LAT+LON
sigl <- ~ sin(LAT/2*pi)+LON
shl <- ~ 1
shl <- ~ LAT+LON
library(ismev)
gev_results <- data_long %>%
  group_by(STATIONS, LAT, LON) %>%
  summarise(gev_fit = list(gev.fit(Precipitation, data_long,
                                   mul = c(2,3), sigl = c(2,3), type = "GEV")))

#######################################################-

### a) Modèle 1 --------
loc.form0 <- ~ coord[,1] + coord[,2]
scale.form0 <- ~ coord[,1] + coord[,2]
shape.form0 <- ~ 1

m1 = fitspatgev(D_Mat, coord, loc.form0, scale.form0, shape.form0, corr =TRUE)
#m1 <- fitspatgev(t(PL[,-c(1,2)]), coord, loc.form0, scale.form0, shape.form0,corr =TRUE)

m1_loc_intercept = m1$param[1]
m1_loc_lat = m1$param[2]
m1_loc_lon = m1$param[3]

m1_scale_intercept = m1$param[4]
m1_scale_lat = m1$param[5]
m1_scale_lon = m1$param[6]

m1_shape_intercept = m1$param[7]


### b) Modèle 2 --------
loc.form1 <- ~ coord[,1] + coord[,2]
scale.form1 <- ~ coord[,1] 
shape.form1 <- ~ 1

m2 = fitspatgev(D_Mat, coord, loc.form1, scale.form1, shape.form1, corr =TRUE)
#m2 <- fitspatgev(t(PL[,-c(1,2)]), coord, loc.form1, scale.form1, shape.form1)

m2_loc_intercept = m2$param[1]
m2_loc_lat = m2$param[2]
m2_loc_lon = m2$param[3]

m2_scale_intercept = m2$param[4]
m2_scale_lat = m2$param[5]

m2_shape_intercept = m2$param[6]



## 2) Modélisation spatio-temporelle avec dépendance sinusoïdale-----------

#years <- rep(1:24, each = ncol(D_Mat))
#temp.cov = matrix(years, 24, byrow=T)
#colnames(temp.cov) =c(1:1568)

years <- rep(1:24)
temp.cov <- data.frame(
  time = years,                     # Temps linéaire (pour loc et scale)
  sin_time = sin(2 * pi * years/5),   # Terme sinusoïdal 
  cos_time = cos(2 * pi * years/5)    # Terme cosinusoïdal
)
### a) Modèle 1------------
loc.form0 <- y ~ lat + lon
scale.form0 <- y ~ lat + lon
shape.form0 <- y ~ 1
temp.form.loc0 <- ~ time
temp.form.scale0 <- ~ time
temp.form.shape0 <- ~ time 

M1 <- fitspatgev(D_Mat, coord, loc.form0, scale.form0, shape.form0, temp.form.loc = temp.form.loc0,
                 temp.form.scale = temp.form.scale0, temp.form.shape = temp.form.shape0,
                 temp.cov = temp.cov, corr =TRUE)

M1_loc_intercept = M1$param[1]
M1_loc_lat = M1$param[2]
M1_loc_lon = M1$param[3]
M1_scale_intercept = M1$param[4]
M1_scale_lat = M1$param[5]
M1_scale_lon = M1$param[6]
M1_shape_intercept = M1$param[7]
M1_temp_loc = M1$param[8]
M1_temp_scale = M1$param[9]
M1_temp_shape = M1$param[10]

### b) Modèle 2------------------
loc.form1 <- y ~ lat + lon
scale.form1 <- y ~ lat 
shape.form1 <- y ~ 1
temp.form.loc1 <- ~ time
temp.form.scale1 <- ~ time
temp.form.shape1 <- ~ time

M2 <- fitspatgev(D_Mat, coord, loc.form1, scale.form1, shape.form1,  temp.form.loc = temp.form.loc1,
                 temp.form.scale = temp.form.scale1, temp.form.shape = temp.form.shape1,
                 temp.cov = temp.cov, corr =TRUE)

M2_loc_intercept = M2$param[1]
M2_loc_lat = M2$param[2]
M2_loc_lon = M2$param[3]
M2_scale_intercept = M2$param[4]
M2_scale_lat = M2$param[5]
M2_shape_intercept = M2$param[6]
M2_temp_loc = M2$param[7]
M2_temp_scale = M2$param[8]
M2_temp_shape = M2$param[9]


### c) Modèle 3------------------
loc.form2 <- y ~ lat + lon
scale.form2 <- y ~ lat + lon
shape.form2 <- y ~ 1
temp.form.loc2 <- ~ time
temp.form.scale2 <- ~ time
temp.form.shape2 <- ~ sin_time + cos_time

M3 <- fitspatgev(D_Mat, coord, loc.form2, scale.form2, shape.form2,  temp.form.loc = temp.form.loc2,
                 temp.form.scale = temp.form.scale2, temp.form.shape = temp.form.shape2,
                 temp.cov = temp.cov, corr =TRUE)

M3_loc_intercept = M3$param[1]
M3_loc_lat = M3$param[2]
M3_loc_lon = M3$param[3]
M3_scale_intercept = M3$param[4]
M3_scale_lat = M3$param[5]
M3_scale_lon = M3$param[6]
M3_shape_intercept = M3$param[7]
M3_temp_loc = M3$param[8]
M3_temp_scale = M3$param[9]
M3_temp_shape_sin = M3$param[10]
M3_temp_shape_cos = M3$param[11]

### d) Modèle 4 ------------------
loc.form3 <- y ~ lat + lon
scale.form3 <- y ~ lat 
shape.form3 <- y ~ 1
temp.form.loc3 <- ~ time
temp.form.scale3 <- ~ time
temp.form.shape3 <- ~ sin_time + cos_time

M4 <- fitspatgev(D_Mat, coord, loc.form3, scale.form3, shape.form3,  temp.form.loc = temp.form.loc3,
                 temp.form.scale = temp.form.scale3, temp.form.shape = temp.form.shape3,
                 temp.cov = temp.cov, corr =TRUE)

M4_loc_intercept = M4$param[1]
M4_loc_lat = M4$param[2]
M4_loc_lon = M4$param[3]
M4_scale_intercept = M4$param[4]
M4_scale_lat = M4$param[5]
M4_shape_intercept = M4$param[6]
M4_temp_loc = M4$param[7]
M4_temp_scale = M4$param[8]
M4_temp_shape_sin = M4$param[9]
M4_temp_shape_cos = M4$param[10]

0### e) Modèle 5 ------------------
loc.form4 <- y ~ lat + lon
scale.form4 <- y ~ lat 
shape.form4 <- y ~ 1
temp.form.loc4 <- ~ time
temp.form.scale4 <- ~ time
temp.form.shape4 <- ~ sin_time 

M5 <- fitspatgev(D_Mat, coord, loc.form4, scale.form4, shape.form4,  temp.form.loc = temp.form.loc4,
                 temp.form.scale = temp.form.scale4, temp.form.shape = temp.form.shape4,
                 temp.cov = temp.cov, corr =TRUE)

M5_loc_intercept = M5$param[1]
M5_loc_lat = M5$param[2]
M5_loc_lon = M5$param[3]
M5_scale_intercept = M5$param[4]
M5_scale_lat = M5$param[5]
M5_shape_intercept = M5$param[6]
M5_temp_loc = M5$param[7]
M5_temp_scale = M5$param[8]
M5_temp_shape_sin = M5$param[9]

### f) Modèle 6------------------
loc.form5 <- y ~ lat + lon
scale.form5 <- y ~ lat + lon
shape.form5 <- y ~ 1
temp.form.loc5 <- ~ time
temp.form.scale5 <- ~ time
temp.form.shape5 <- ~ sin_time 

M6 <- fitspatgev(D_Mat, coord, loc.form5, scale.form5, shape.form5,  temp.form.loc = temp.form.loc5,
                 temp.form.scale = temp.form.scale5, temp.form.shape = temp.form.shape5,
                 temp.cov = temp.cov, corr =TRUE)

M6_loc_intercept = M6$param[1]
M6_loc_lat = M6$param[2]
M6_loc_lon = M6$param[3]
M6_scale_intercept = M6$param[4]
M6_scale_lat = M6$param[5]
M6_scale_lon = M6$param[6]
M6_shape_intercept = M6$param[7]
M6_temp_loc = M6$param[8]
M6_temp_scale = M6$param[9]
M6_temp_shape_sin = M6$param[10]

### g) Modèle 7 ------------------
loc.form6 <- y ~ lat + lon
scale.form6 <- y ~ lat 
shape.form6 <- y ~ 1
temp.form.loc6 <- ~ time
temp.form.scale6 <- ~ time
temp.form.shape6 <- ~ cos_time 

M7 <- fitspatgev(D_Mat, coord, loc.form6, scale.form6, shape.form6,  temp.form.loc = temp.form.loc6,
                 temp.form.scale = temp.form.scale6, temp.form.shape = temp.form.shape6,
                 temp.cov = temp.cov, corr =TRUE)

M7_loc_intercept = M7$param[1]
M7_loc_lat = M7$param[2]
M7_loc_lon = M7$param[3]
M7_scale_intercept = M7$param[4]
M7_scale_lat = M7$param[5]
M7_shape_intercept = M7$param[6]
M7_temp_loc = M7$param[7]
M7_temp_scale = M7$param[8]
M7_temp_shape_cos = M7$param[9]

### h) Modèle 8------------------
loc.form7 <- y ~ lat + lon
scale.form7 <- y ~ lat + lon
shape.form7 <- y ~ 1
temp.form.loc7 <- ~ time
temp.form.scale7 <- ~ time
temp.form.shape7 <- ~ cos_time 

M8 <- fitspatgev(D_Mat, coord, loc.form7, scale.form7, shape.form7,  temp.form.loc = temp.form.loc7,
                 temp.form.scale = temp.form.scale7, temp.form.shape = temp.form.shape7,
                 temp.cov = temp.cov, corr =TRUE)

M8_loc_intercept = M8$param[1]
M8_loc_lat = M8$param[2]
M8_loc_lon = M8$param[3]
M8_scale_intercept = M8$param[4]
M8_scale_lat = M8$param[5]
M8_scale_lon = M8$param[6]
M8_shape_intercept = M8$param[7]
M8_temp_loc = M8$param[8]
M8_temp_scale = M8$param[9]
M8_temp_shape_cos = M8$param[10]


## 3) Modélisation spatio-temporelle sans dépendance sinusoïdale-----------
### a- Possibilité que xi soit constant ds le temps ----------


#years <- rep(1:24, each = ncol(D_Mat))
#temp.cov = matrix(years, 24, byrow=T)
#colnames(temp.cov) =c(1:1568)

years <- rep(1:24)
temp.cov2 <- data.frame(
  time = years    # Terme cosinusoïdal
)
### a) Modèle 1------------
loc.form_0 <- y ~ lat + lon
scale.form_0 <- y ~ lat + lon
shape.form_0 <- y ~ 1
temp.form.loc_0 <- ~ time
temp.form.scale_0 <- ~ time
#temp.form.shape_0 <- ~ rep(mean(PSI_Temps), length(PSI_Temps)) 

Mod1 <- fitspatgev(D_Mat, coord, loc.form_0, scale.form_0, shape.form_0, temp.form.loc = temp.form.loc_0,
                   temp.form.scale = temp.form.scale_0, #temp.form.shape = temp.form.shape_0,
                   temp.cov = temp.cov2, corr =TRUE)

Mod1_loc_intercept = Mod1$param[1]
Mod1_loc_lat = Mod1$param[2]
Mod1_loc_lon = Mod1$param[3]
Mod1_scale_intercept = Mod1$param[4]
Mod1_scale_lat = Mod1$param[5]
Mod1_scale_lon = Mod1$param[6]
Mod1_shape_intercept = Mod1$param[7]
Mod1_temp_loc = Mod1$param[8]
Mod1_temp_scale = Mod1$param[9]
#Mod1_temp_shape = Mod1$param[10]

### b) Modèle 2------------------
loc.form_1 <- y ~ lat + lon
scale.form_1 <- y ~ lat 
shape.form_1 <- y ~ 1
temp.form.loc_1 <- ~ time
temp.form.scale_1 <- ~ time
#temp.form.shape_1 <- ~ rep(mean(PSI_Temps), length(PSI_Temps)) 

Mod2 <- fitspatgev(D_Mat, coord, loc.form_1, scale.form_1, shape.form_1,  temp.form.loc = temp.form.loc_1,
                   temp.form.scale = temp.form.scale_1, #temp.form.shape = temp.form.shape_1,
                   temp.cov = temp.cov2, corr =TRUE)

Mod2_loc_intercept = Mod2$param[1]
Mod2_loc_lat = Mod2$param[2]
Mod2_loc_lon = Mod2$param[3]
Mod2_scale_intercept = Mod2$param[4]
Mod2_scale_lat = Mod2$param[5]
Mod2_shape_intercept = Mod2$param[6]
Mod2_temp_loc = Mod2$param[7]
Mod2_temp_scale = Mod2$param[8]
#Mod2_temp_shape = Mod2$param[9]


### c) Modèle 3------------------
loc.form_2 <- y ~ lat + lon
scale.form_2 <- y ~ lat + lon
shape.form_2 <- y ~ 1
temp.form.loc_2 <- ~ time
temp.form.scale_2 <- ~ time
temp.form.shape_2 <- ~ time

Mod3 <- fitspatgev(D_Mat, coord, loc.form_2, scale.form_2, shape.form_2,  temp.form.loc = temp.form.loc_2,
                   temp.form.scale = temp.form.scale_2, temp.form.shape = temp.form.shape_2,
                   temp.cov = temp.cov2, corr =TRUE)

Mod3_loc_intercept = Mod3$param[1]
Mod3_loc_lat = Mod3$param[2]
Mod3_loc_lon = Mod3$param[3]
Mod3_scale_intercept = Mod3$param[4]
Mod3_scale_lat = Mod3$param[5]
Mod3_scale_lon = Mod3$param[6]
Mod3_shape_intercept = Mod3$param[7]
Mod3_temp_loc = Mod3$param[8]
Mod3_temp_scale = Mod3$param[9]
Mod3_temp_shape = Mod3$param[10]


### d) Modèle 4 ------------------
loc.form_3 <- y ~ lat + lon
scale.form_3 <- y ~ lat 
shape.form_3 <- y ~ 1
temp.form.loc_3 <- ~ time
temp.form.scale_3 <- ~ time
temp.form.shape_3 <- ~ time

Mod4 <- fitspatgev(D_Mat, coord, loc.form_3, scale.form_3, shape.form_3,  temp.form.loc = temp.form.loc_3,
                   temp.form.scale = temp.form.scale_3, temp.form.shape = temp.form.shape_3,
                   temp.cov = temp.cov2, corr =TRUE)

Mod4_loc_intercept = Mod4$param[1]
Mod4_loc_lat = Mod4$param[2]
Mod4_loc_lon = Mod4$param[3]
Mod4_scale_intercept = Mod4$param[4]
Mod4_scale_lat = Mod4$param[5]
Mod4_shape_intercept = Mod4$param[6]
Mod4_temp_loc = Mod4$param[7]
Mod4_temp_scale = Mod4$param[8]
Mod4_temp_shape = Mod4$param[9]

### a- xi n'est pas constant ds le temps ----------



#years <- rep(1:24, each = ncol(D_Mat))
#temp.cov = matrix(years, 24, byrow=T)
#colnames(temp.cov) =c(1:1568)

years <- rep(1:24)
temp.cov2 <- data.frame(
  time = years    # Terme cosinusoïdal
)

### a) Modèle 1------------------
loc.form_Mod_1 <- y ~ lat + lon
scale.form_Mod_1 <- y ~ lat + lon
shape.form_Mod_1 <- y ~ 1
temp.form.loc_Mod_1 <- ~ time
temp.form.scale_Mod_1 <- ~ time
temp.form.shape_Mod_1 <- ~ time

Mod_1 <- fitspatgev(D_Mat, coord, loc.form_Mod_1, scale.form_Mod_1, shape.form_Mod_1,  temp.form.loc = temp.form.loc_Mod_1,
                    temp.form.scale = temp.form.scale_Mod_1, temp.form.shape = temp.form.shape_Mod_1,
                    temp.cov = temp.cov2, corr =TRUE)

Mod_1_loc_intercept = Mod3$param[1]
Mod_1_loc_lat = Mod3$param[2]
Mod_1_loc_lon = Mod3$param[3]
Mod_1_scale_intercept = Mod3$param[4]
Mod_1_scale_lat = Mod3$param[5]
Mod_1_scale_lon = Mod3$param[6]
Mod_1_shape_intercept = Mod3$param[7]
Mod_1_temp_loc = Mod3$param[8]
Mod_1_temp_scale = Mod3$param[9]
Mod_1_temp_shape = Mod3$param[10]


### b) Modèle 2 ------------------
loc.form_Mod_2 <- y ~ lat + lon
scale.form_Mod_2 <- y ~ lat 
shape.form_Mod_2 <- y ~ 1
temp.form.loc_Mod_2 <- ~ time
temp.form.scale_Mod_2 <- ~ time
temp.form.shape_Mod_2 <- ~ time

Mod_2 <- fitspatgev(D_Mat, coord, loc.form_Mod_2, scale.form_Mod_2, shape.form_Mod_2,  temp.form.loc = temp.form.loc_Mod_2,
                    temp.form.scale = temp.form.scale_Mod_2, temp.form.shape = temp.form.shape_Mod_2,
                    temp.cov = temp.cov2, corr =TRUE)

Mod_2_loc_intercept = Mod4$param[1]
Mod_2_loc_lat = Mod4$param[2]
Mod_2_loc_lon = Mod4$param[3]
Mod_2_scale_intercept = Mod4$param[4]
Mod_2_scale_lat = Mod4$param[5]
Mod_2_shape_intercept = Mod4$param[6]
Mod_2_temp_loc = Mod4$param[7]
Mod_2_temp_scale = Mod4$param[8]
Mod_2_temp_shape = Mod4$param[9]

# C- CHOIX DU MEILLEUR MODELE----------
# Modèles avec dépenance sinusoïdale ----------
## Deviance ------------
dev1 = m1$deviance
dev2 = m2$deviance

dev3 = M1$deviance
dev4 = M2$deviance
dev5 = M3$deviance
dev6 = M4$deviance
dev7 = M5$deviance
dev8 = M6$deviance
dev9 = M7$deviance
dev10 = M8$deviance

tab_deviance = data.frame(Modeles = c("m1", "m2", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"),
                          Deviances = c(dev1, dev2, dev3, dev4, dev5, dev6, dev7, dev8, dev9, dev10))
tab_deviance

#Les modèles 6 et 7 ont les plus faibles déviances.
## AIC ----------
aic_funct = function(modele){
  k = length(modele$param)
  aic_mod = -2*modele$logLik + 2*k
  cat("\nlogLik =", modele$logLik, "; k = ", k)
  cat("\nAIC = ", aic_mod)
  return(aic_mod)
}
AIC1 = aic_funct(m1)
AIC2 = aic_funct(m2)
AIC3 = aic_funct(M1)
AIC4 = aic_funct(M2)
AIC5 = aic_funct(M3)
AIC6 = aic_funct(M4)
AIC7 = aic_funct(M5)
AIC8 = aic_funct(M6)
AIC9 = aic_funct(M7)
AIC10 = aic_funct(M8)


## BIC ---------

bic_funct <- function(modele, n) {
  k <- length(modele$param)
  bic_mod <- -2 * modele$logLik + k * log(n)  # n = nombre d'observations
  return(bic_mod)
}

n = 24 * ncol(D_Mat)
BIC1 = bic_funct(m1, n)
BIC2 = bic_funct(m2, n)
BIC3 = bic_funct(M1, n)
BIC4 = bic_funct(M2, n)
BIC5 = bic_funct(M3, n)
BIC6 = bic_funct(M4, n)
BIC7 = bic_funct(M5, n)
BIC8 = bic_funct(M6, n)
BIC9 = bic_funct(M7, n)
BIC10 = bic_funct(M8, n)


tab_critere = data.frame(Modeles = c("m1", "m2", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8"),
                         Deviances = c(dev1, dev2, dev3, dev4, dev5, dev6, dev7, dev8, dev9, dev10),
                         AIC = c(AIC1, AIC2, AIC3, AIC4, AIC5, AIC6, AIC7, AIC8, AIC9, AIC10),
                         BIC = c(BIC1, BIC2, BIC3, BIC4, BIC5, BIC6, BIC7, BIC8, BIC9, BIC10))
tab_critere

# Modèles sans dépenance sinusoïdale ----------
## Deviance ------------
dev1 = m1$deviance
dev2 = m2$deviance
devMod1 = Mod1$deviance
devMod2 = Mod2$deviance
devMod3 = Mod3$deviance
devMod4 = Mod4$deviance

AIC1 = aic_funct(m1)
AIC2 = aic_funct(m2)
AICMod1 = aic_funct(Mod1)
AICMod2 = aic_funct(Mod2)
AICMod3 = aic_funct(Mod3)
AICMod4 = aic_funct(Mod4)

n = 24 * ncol(D_Mat)
BIC1 = bic_funct(m1, n)
BIC2 = bic_funct(m2, n)
BICMod1 = bic_funct(Mod1, n)
BICMod2 = bic_funct(Mod2, n)
BICMod3 = bic_funct(Mod3, n)
BICMod4 = bic_funct(Mod4, n)
tab_critere2 = data.frame(Modeles = c("m1", "m2", "Mod1", "Mod2", "Mod3", "Mod4"),
                          Deviances = c(dev1, dev2, devMod1, devMod2, devMod3, devMod4),
                          AIC = c(AIC1, AIC2, AICMod1, AICMod2, AICMod3, AICMod4),
                          BIC = c(BIC1, BIC2, BICMod1, BICMod2, BICMod3, BICMod4))
tab_critere2

# Modèles avec xi lin et sans dépendance sinusoïdale ----------
## Deviance ------------
dev1 = m1$deviance
dev2 = m2$deviance
devMod_1 = Mod_1$deviance
devMod_2 = Mod_2$deviance


AIC1 = aic_funct(m1)
AIC2 = aic_funct(m2)
AICMod_1 = aic_funct(Mod_1)
AICMod_2 = aic_funct(Mod_2)


n = 24 * ncol(D_Mat)
BIC1 = bic_funct(m1, n)
BIC2 = bic_funct(m2, n)
BICMod_1 = bic_funct(Mod_1, n)
BICMod_2 = bic_funct(Mod_2, n)

tab_critere3 = data.frame(Modeles = c("m1", "m2", "Mod_1", "Mod_2"),
                          Deviances = c(dev1, dev2, devMod_1, devMod_2),
                          AIC = c(AIC1, AIC2, AICMod_1, AICMod_2),
                          BIC = c(BIC1, BIC2, BICMod_1, BICMod_2))
tab_critere3

critères = list(tab_critere, tab_critere2, tab_critere3)
# D- DIAGNOSTIQUE DU MEILLEUR MODELE--------
## 0) qqgev --------

qqgev(M1)
qqgev(m1)
qqgev(Mod1)


# E- NIVEAUX DE RETOUR ----------
## 1) En fonction de lat, lon, temps et période de retour--------
Niveau_Retour= function(la, lo, t, Annee){
  time=0 ; z_T = 0
  if (Annee <=2023 & Annee >=2000){
    time =(Annee-2000+1)
    
    #Meilleur modèle spatial
    #mu = m1_loc_intercept + m1_loc_lat*la + m1_loc_lon*lo 
    #sigma = m1_scale_intercept + m1_scale_lat*la + m1_scale_lon*lo 
    #xi = m1_shape_intercept 
    
    
    #Meilleur modèle spatio temporel sans dépendance sinusoïdale et xi dpt uniqmt du tps
    #mu = Mod_1_loc_intercept + Mod_1_loc_lat*la + Mod_1_loc_lon*lo + Mod_1_temp_loc*time 
    #sigma = Mod_1_scale_intercept + Mod_1_scale_lat*la + Mod_1_scale_lon*lo + Mod_1_temp_scale*time
    #xi = Mod_1_shape_intercept
    
    #mu = Mod3_loc_intercept + Mod3_loc_lat*la + Mod3_loc_lon*lo + Mod3_temp_loc*time 
    #sigma = Mod3_scale_intercept + Mod3_scale_lat*la + Mod3_scale_lon*lo + Mod3_temp_scale*time
    #xi = Mod3_shape_intercept
    
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
    mu = M1_loc_intercept + M1_loc_lat*la + M1_loc_lon*lo + M1_temp_loc*time 
    sigma = M1_scale_intercept + M1_scale_lat*la + M1_scale_lon*lo + M1_temp_scale*time
    xi = M1_shape_intercept +  M1_temp_shape*time
    
    
    #Meilleur modèle spatio temporel avec normalisation de lat et lon et periode 24
    #mu = M_sc_loc_intercept + M_sc_loc_lat*la + M_sc_loc_lon*lo + M_sc_temp_loc*time 
    #sigma = M_sc_scale_intercept + M_sc_scale_lat*la + M_sc_scale_lon*lo + M_sc_temp_scale*time
    #xi = M_sc_shape_intercept +  M_sc_temp_shape_sin*time
    
    if (abs(xi) <= 1e-6) {
      # Cas Gumbel : xi = 0
      z_T <- mu - sigma * log(-log(1 - 1 / t))
    } else {
      # Cas général : xi != 0
      z_T <- mu + (sigma/xi) * ( (-log(1 - 1 / t))^(-xi) - 1)
    }
  }
  
  list(MU=mu, Sig= sigma, Rtn =z_T)
}

## 2) En fonction de la période de retour -------
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



#map(M6, x = seq(-16.75, 10.75, 0.5), y = seq(4.25, 17.75, 0.5), covariates = NULL,
#    param = "quant", ret.per = 100, col = terrain.colors(64), plot.contour = TRUE)

save.image("analyse_exploratoire.rds")
