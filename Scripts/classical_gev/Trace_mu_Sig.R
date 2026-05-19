load("trace_mu_sig.rds")

library(SpatialExtremes); library(RColorBrewer) ; library(ggplot2); library(latex2exp); library(patchwork)


# VISUALISATION ANALYSE EXPLORATOIRE -----
#1 Définition des bases de données des paramètres estimés. -----
## a Data latitude --------
Bande_Latitude  <- data.frame(
  Latitude = rep(Bande_Lat,24),
  Mu = as.vector(MU_Lat), 
  Sig = as.vector(SIG_Lat),
  Year  = sort(c(rep(2000,length(Bande_Lat)), rep(2001,length(Bande_Lat)),  rep(2002,length(Bande_Lat)), rep(2003,length(Bande_Lat)),  
rep(2004,length(Bande_Lat)), rep(2005,length(Bande_Lat)),  rep(2006,length(Bande_Lat)), rep(2007,length(Bande_Lat)),  
rep(2008,length(Bande_Lat)), rep(2009,length(Bande_Lat)),  rep(2010,length(Bande_Lat)), rep(2011,length(Bande_Lat)),  
rep(2012,length(Bande_Lat)), rep(2013,length(Bande_Lat)),  rep(2014,length(Bande_Lat)), rep(2015,length(Bande_Lat)),    
rep(2016,length(Bande_Lat)), rep(2017,length(Bande_Lat)),  rep(2018,length(Bande_Lat)), rep(2019,length(Bande_Lat)),  
rep(2020,length(Bande_Lat)), rep(2021,length(Bande_Lat)),  rep(2022,length(Bande_Lat)), rep(2023,length(Bande_Lat))),T) 
)

#View(Bande_Latitude)

PSI_Latitude= PSI_Lat; Annee= 2023:2000

## b Data longitude --------
Bande_Longitude  <- data.frame(
  Longitude  = rep(Bande_Lon,24),
  Mu = as.vector(MU_Lon), 
  Sig = as.vector(SIG_Lon),
  Year  = sort(c(rep(2000,length(Bande_Lon)), rep(2001,length(Bande_Lon)),  rep(2002,length(Bande_Lon)), rep(2003,length(Bande_Lon)),  
rep(2004,length(Bande_Lon)), rep(2005,length(Bande_Lon)),  rep(2006,length(Bande_Lon)), rep(2007,length(Bande_Lon)),  
rep(2008,length(Bande_Lon)), rep(2009,length(Bande_Lon)),  rep(2010,length(Bande_Lon)), rep(2011,length(Bande_Lon)),  
rep(2012,length(Bande_Lon)), rep(2013,length(Bande_Lon)),  rep(2014,length(Bande_Lon)), rep(2015,length(Bande_Lon)),    
rep(2016,length(Bande_Lon)), rep(2017,length(Bande_Lon)),  rep(2018,length(Bande_Lon)), rep(2019,length(Bande_Lon)),  
rep(2020,length(Bande_Lon)), rep(2021,length(Bande_Lon)),  rep(2022,length(Bande_Lon)), rep(2023,length(Bande_Lon))),T) 
)

#View(Bande_Longitude)

PSI_Longitude  = PSI_Lon

## c Data temps -----------

Bande_Temps  <- data.frame(
  Year = seq(2023, 2000, -1),
  Mu = MU_Temps, 
  Sig = SIG_Temps,
  Psi = PSI_Temps
)
#View(Bande_Temps)
#2 Visualisation des graphiques-----------


## a Paramètres gev en fonction de la latitude ------------
ggplot(Bande_Latitude , aes(x = Latitude, y = Mu)) +
 geom_point(colour = "blue", size = 2) +  
facet_wrap(~ Year) +

scale_x_continuous(name = "Latitude", 
                     breaks = seq(4, 18, by = 2),
                     labels = function(x) paste0(x, "°N"))+

theme(legend.title = element_text(hjust = 0.5, size = 10, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 12, ),
  strip.text = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 10, face = "bold"),
  axis.text.y = element_text(size = 10, face = "bold"),
  axis.title.x = element_text(size = 12, face = "bold"), 
  axis.title.y = element_text(size = 12, face = "bold"))+

  labs(title = TeX("$mu_{lat}$"), y = TeX("$mu$") )




ggplot(Bande_Latitude , aes(x = Latitude, y = Sig)) +
 geom_point( colour = "red", size = 2) +  
facet_wrap(~ Year) +

scale_x_continuous(name = "Latitude", 
                     breaks = seq(8, 14, by = 2),
                     labels = function(x) paste0(x, "°N"))+

theme(legend.title = element_text(hjust = 0.5, size = 10, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 12, ),
  strip.text = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 10, face = "bold"),
  axis.text.y = element_text(size = 10, face = "bold"),
  axis.title.x = element_text(size = 12, face = "bold"), 
  axis.title.y = element_text(size = 12, face = "bold"))+

  labs(title = TeX("$sigma_{lat}$"),
  y = TeX("$sigma$"))



## b Paramètres gev en fonction de la longitude ------------

ggplot(Bande_Longitude , aes(x = Longitude , y = Mu)) +
 geom_point( colour = "blue", size = 2) +  
facet_wrap(~ Year) +

scale_x_continuous(name = "Longitude", 
                   breaks = seq(-17, 11, by = 4),
                   labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E")))+

theme(legend.title = element_text(hjust = 0.5, size = 10, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 12),
  strip.text = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 10, face = "bold"),
  axis.text.y = element_text(size = 10, face = "bold"),
  axis.title.x = element_text(size = 12, face = "bold"), 
  axis.title.y = element_text(size = 12, face = "bold"))+

  labs(title = TeX("$mu_{lon}$"),
  y = TeX("$mu$"))



ggplot(Bande_Longitude , aes(x = Longitude, y = Sig)) +
 geom_point( colour = "red", size = 2) +  
facet_wrap(~ Year) +

scale_x_continuous(name = "Longitude", 
                   breaks = seq(-17, 11, by = 4),
                   labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E")))+

theme(legend.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 12, ),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 10, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 12, face = "bold"), 
      axis.title.y = element_text(size = 12, face = "bold"))+

  labs(title = TeX("$sigma_{lon}$"),
  y = TeX("$sigma$"))

## c Paramètres gev en fonction du temps ------------
fig1 = ggplot(Bande_Temps , aes(x = Year, y = Sig)) +
  geom_line(colour = "red", size = 2) +  
  
  scale_x_continuous(name = "Years")+
  
  theme(legend.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, ),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"))+
  
  labs(title = TeX("$sigma_t$"),
       y = TeX("$sigma$"))


fig2 = ggplot(Bande_Temps , aes(x = Year, y = Mu)) +
  geom_line(colour = "red", size = 2) +  
  
  scale_x_continuous(name = "Years")+
  
  theme(legend.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, ),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"))+
  
  labs(title = TeX("$mu_t$"),
       y = TeX("$mu$"))


fig3 = ggplot(Bande_Temps , aes(x = Year, y = Psi)) +
  geom_line(colour = "red", size = 2) +  
  
  scale_x_continuous(name = "Years")+
  
  theme(legend.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        plot.title = element_text(hjust = 0.5, size = 12, ),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"))+
  
  labs(title = TeX("$xi_t$"),
       y = TeX("$xi$"))

fig1 / fig2 / fig3

#VISUALISATION MODELISATION   ---------- 


PL <- read.delim("../AO.txt")
#View(PL)
attach(PL)

coord <- matrix(c(LAT,LON), ncol=2) ;  colnames(coord) <- c("lat", "lon")

D_Mat = matrix(c(A2000, A2001, A2002, A2003, A2004, A2005, A2006, A2007, A2008, A2009,
            A2010, A2011, A2012, A2013, A2014 ,A2015 ,A2016, A2017, A2018, A2019, A2020, A2021, A2022 , A2023),
              ncol=1568, byrow=TRUE )

#View(D_Mat)



RESULTAT=Result(50) # 50-year RL

Variation_Mu_Sig  <- data.frame(
  Latitude = RESULTAT$Mu_mat[,1],
  Longitude = RESULTAT$Mu_mat[,2],
  Mu = RESULTAT$Mu_mat[,3], 
  Sig = RESULTAT$Sig_mat[,3],
  Year  =  RESULTAT$Mu_mat[,4],
Return_level = RESULTAT$Return.level[,3]
 )
attach(Variation_Mu_Sig)
#View(Variation_Mu_Sig)


ggplot() +
  # Couleur des frontières
  geom_point(data = Variation_Mu_Sig, aes(x = Longitude, y = Latitude, color =  Return_level), size = 8, stat = "identity") + 
  #scale_color_gradient(low = c("darkblue","lightblue"), high = c("red")) + # Dégradé de couleur
  
  scale_color_gradient(low = c("lightblue","darkblue", "green"), high = c("yellow", "orange", "red")) + # Dégradé de couleur
  #scale_fill_viridis_d(option = "A")+
  
  scale_x_continuous(name = "Longitude", 
                     breaks = seq(-17, 11, by = 6),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(name = "Latitude", 
                     breaks = seq(4, 20, by = 3),
                     labels = function(y) paste0(y, "°N")) +
  
  labs(title = "50-years return level", 
       #  x = "Longitude",
       #    y = "Latitude",
       color = "(mm/jr)")+
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

## Les précipitations selon les classes

Variation_Mu_Sig$Return_level_cut = cut(Return_level,
                       breaks = 6,
                       labels = c("N1","N2","N3", "N4","N5","N6") )

#table(Return_level_cat)

##Niveau de retour sur ########
ggplot() +
 # Couleur des frontières
  geom_point(data = Variation_Mu_Sig, aes(x = Longitude, y = Latitude, color =  Return_level_cut), size =5,stat = "identity") + 
 #scale_color_gradient(low = c("darkblue","lightblue"), high = c("red")) + # Dégradé de couleur
  scale_color_manual(values = c("N1" = "darkblue", 
                                "N2" = "lightblue",
                                "N3" = "green",
                                "N4" = "yellow",
                                "N5" = "orange",
                                "N6" = "red")) +
 # scale_color_gradient(low = c("darkblue","lightblue", "green"), high = c("yellow", "orange", "red")) + # Dégradé de couleur
#scale_fill_viridis_d(option = "A")+

  scale_x_continuous(name = "Longitude", 
                     breaks = seq(-17, 11, by = 6),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(name = "Latitude", 
                     breaks = seq(4, 20, by = 3),
                     labels = function(y) paste0(y, "°N")) +

 labs(title = "10-year return level", 
     #  x = "Longitude",
   #    y = "Latitude",
       color = "Level(mm/jr)")+
  theme_minimal() +
  facet_wrap(~ Year) +

borders("world", regions = c("Benin", "Burkina Faso", "Ivory-Coast", 
               "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
                "Liberia", "Mali", "Mauritania", "Niger", 
                "Nigeria", "Senegal", "Sierra Leone", "Togo"), 
         colour = "black",size=.8) +

theme(legend.title = element_text(hjust = 0.5, size = 10, face = "bold", ), 
      legend.position="top", legend.text = element_text(size = 10, face = "bold"),
#legend.position.inside =c(15,18), 
plot.title = element_text(hjust = 0.5, size = 12, ), 
strip.text = element_text(size = 8, face = "bold"), 
axis.text.x = element_text(size = 10, face = "bold"),
axis.text.y = element_text(size = 10, face = "bold"),
axis.title.x = element_text(size = 12, face = "bold"), 
axis.title.y = element_text(size = 12, face = "bold"))+
  
  coord_fixed( xlim = c(-16.75, 10.1), ylim = c(4.25, 17.75))  

## Niveau de retour sur 10 et 100 ans pr les 10 premieres et dernieres annees########


RESULTAT_50=Result(t=50) ; RESULTAT_100=Result(t=100)
 

Premier_10  <- data.frame(
  Latitude = c(RESULTAT_50$Mu_mat[1:15680,1] ,RESULTAT_100$Mu_mat[1:15680,1] ),
  Longitude = c(RESULTAT_50$Mu_mat[1:15680,2] ,RESULTAT_100$Mu_mat[1:15680,2] ),
  Year  =  c(RESULTAT_50$Mu_mat[1:15680,4], RESULTAT_100$Mu_mat[1:15680,4] ),
Return_level = c(RESULTAT_50$Return.level[1:15680,3], RESULTAT_100$Return.level[1:15680,3]),
 periode=c(rep(" 50-years",15680),rep("100-years",15680))
 )


Dernier_10  <- data.frame(
  Latitude = c(RESULTAT_50$Mu_mat[-c(1:21952),1] ,RESULTAT_100$Mu_mat[-c(1:21952),1] ),
  Longitude = c(RESULTAT_50$Mu_mat[-c(1:21952),2] ,RESULTAT_100$Mu_mat[-c(1:21952),2] ),
  Year  =  c(RESULTAT_50$Mu_mat[-c(1:21952),4], RESULTAT_100$Mu_mat[-c(1:21952),4] ),
Return_level = c(RESULTAT_50$Return.level[-c(1:21952),3], RESULTAT_100$Return.level[-c(1:21952),3]),
 periode=c(rep(" 50-years",15680),rep("100-years",15680))
 )



ggplot() +
 # Couleur des frontières
  geom_point(data = Premier_10, aes(x = Longitude, y = Latitude, color =  Return_level), size =8,stat = "identity") + # Points colorés par température
 #scale_color_gradient(low = c("darkblue","lightblue"), high = c("red")) + # Dégradé de couleur

scale_color_gradient(low = c("lightblue","darkblue", "green"), high = c("yellow", "orange", "red")) + # Dégradé de couleur
#scale_fill_viridis_d(option = "A")+

  scale_x_continuous(name = "Longitude", 
                     breaks = seq(-16, 10, by = 8),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(name = "Latitude", 
                     breaks = seq(4, 20, by = 3),
                     labels = function(y) paste0(y, "°N")) +

 labs(title = "Comparison of the 50-years and 100-years return level values", 
     #  x = "Longitude",
   #    y = "Latitude",
       color = "Value(mm/jr)")+
  theme_minimal() +
  #facet_wrap(~ Year) +
facet_grid(rows = vars(periode), cols = vars(Year)) 

#borders("world", regions = c("Benin", "Burkina Faso", "Côte-d'Ivoire", 
 #                "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
 #                "Liberia", "Mali", "Mauritania", "Niger", 
 #                "Nigeria", "Senegal", "Sierra Leone", "Togo"), 
 #         colour = "black",size=.8) +

theme(legend.title = element_text(hjust = 0.5, size = 17, face = "bold"), legend.position="right",
#legend.position.inside =c(15,18), 
plot.subtitle = element_text(hjust = 15, size = 30, face = "bold"),
plot.title = element_text(hjust = 0.5, size = 18, ),
strip.text = element_text(size = 14, face = "bold"),
axis.text.x = element_text(size = 12, face = "bold"),
axis.text.y = element_text(size = 12, face = "bold"),
axis.title.x = element_text(size = 17, face = "bold"), 
axis.title.y = element_text(size = 17, face = "bold"))+
  coord_fixed( xlim = c(-16.75, 10.1), ylim = c(4.25, 17.75))  



ggplot() +
  # Couleur des frontières
  geom_point(data = Dernier_10, aes(x = Longitude, y = Latitude, color =  Return_level), size =8,stat = "identity") + # Points colorés par température
  #scale_color_gradient(low = c("darkblue","lightblue"), high = c("red")) + # Dégradé de couleur
  
  scale_color_gradient(low = c("lightblue","darkblue", "green"), high = c("yellow", "orange", "red")) + # Dégradé de couleur
  #scale_fill_viridis_d(option = "A")+
  
  scale_x_continuous(name = "Longitude", 
                     breaks = seq(-16, 10, by = 8),
                     labels = function(x) paste0(abs(x), ifelse(x < 0, "°W", "°E"))) + # Ajout de W ou E
  scale_y_continuous(name = "Latitude", 
                     breaks = seq(4, 20, by = 3),
                     labels = function(y) paste0(y, "°N")) +
  
  labs(title = "Comparison of the 50-years and 100-years return level values", 
       #  x = "Longitude",
       #    y = "Latitude",
       color = "Value(mm/jr)")+
  theme_minimal() +
  #facet_wrap(~ Year) +
  facet_grid(rows = vars(periode), cols = vars(Year))

#borders("world", regions = c("Benin", "Burkina Faso", "Côte-d'Ivoire", 
#                "Gambia", "Ghana", "Guinea", "Guinea-Bissau", 
#                "Liberia", "Mali", "Mauritania", "Niger", 
#                "Nigeria", "Senegal", "Sierra Leone", "Togo"), 
#         colour = "black",size=.8) +

theme(legend.title = element_text(hjust = 0.5, size = 17, face = "bold"), legend.position="right",
      #legend.position.inside =c(15,18), 
      plot.subtitle = element_text(hjust = 15, size = 30, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 18, ),
      strip.text = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 17, face = "bold"), 
      axis.title.y = element_text(size = 17, face = "bold"))+
  coord_fixed( xlim = c(-16.75, 10.1), ylim = c(4.25, 17.75))  



save.image("trace_mu_sig.rds")