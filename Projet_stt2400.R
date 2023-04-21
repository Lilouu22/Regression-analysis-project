#Noms : Lisa Pamela Valente & Harkavaljot Kaur Chhina
# Matricules : 20180163 & 20157163
# Date  de remise : 14 avril 2023
# Projet stt2400
#########################################################################
#Code inspiré des notes de cours de Alejandro (Alexandro) Murua en stt2400

#Import data
library(readr)
maquereau <- read_csv("C:/Users/lisav/OneDrive/Desktop/STT2400/maquereau.csv")
View(maquereau)
head(maquereau)
#Pour utiliser le nom des variables
attach(maquereau)

#On constate qu'il y a aucune donnée catégorielles ou binaires
#Pas de transformations à faire de ce côté
class(maquereau)

#combien de modèles possibles? 2^10
#Données sont dans des unités qui concordent
# Y cest egg. densité

egg <- subset(maquereau,select=-c(...1))
pairs(egg)
cor(egg)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}


pairs(egg, panel = panel.smooth, pch=16, upper.panel = NULL, diag.panel = panel.hist)
#modèle additif
m0<-lm(egg.dens ~ b.depth+lat+lon+time+flow+s.depth+temp.surf+temp.20m+net.area+c.dist, data = egg)
summary(m0)

#On va télécharger toutes les librairies nécessaires 
library(olsrr)
library(MASS)
library(lmtest)
library(zoo)
library(stats)
library(performance)
library(caTools)
library(car)
library(quantmod)
library(corrplot)
vif(m0)

#VIF = 15,427 109  temp.surf
m1<-lm(egg.dens ~ b.depth+lat+lon+time+flow+s.depth+temp.20m+net.area+c.dist, data = egg)
summary(m1)
vif(m1) # ON ENLÈVE LON = 8,044456
m2<-lm(egg.dens ~ b.depth+lat+time+flow+s.depth+temp.20m+net.area+c.dist, data = egg)
summary(m2)
vif(m2)  #VIF  temp.20m = 5,969322 car on a pris comme décision que cétait too much
#plot(m2)
egg_start<-lm(egg.dens ~ b.depth+lat+time+flow+s.depth+net.area+c.dist, data = egg)
summary(egg_start)
vif(egg_start)  # TOUS EN BAS DE 5
#plot(egg_start) enlever en commentaire pour QQplot, fitted value, scale-location & levrage

#Résidus du modèle additif
residual0<-resid(egg_start)
plot(residual0, ylab = "résidus")
abline(h=0)

###Cook distance 
library(olsrr)
ols_plot_cooksd_bar(egg_start) #312

######Test transformation #########
# Pour la variable réponse 
logYm0 <- mean(log(egg.dens)) 
#logYm0= 3.344241

egg$Gy = egg$egg.dens*( log(egg$egg.dens)-logYm0-1)
Gy.reg = lm(Gy ~ b.depth+lat+time+flow+s.depth+net.area+c.dist, data = egg )

#Graphique de la variable ajoutée
plot(Gy.reg$residuals, egg_start$residuals, xlab = "Modèle avec Gy", ylab="Modèle additif")
Gy.lowess = lowess(egg_start$residuals ~ Gy.reg$residuals , f= 0.75)
lines( Gy.lowess, col="red", lty=2, lwd=2 )

egg.Gy.lm = lm( egg.dens ~ b.depth+lat+time+flow+s.depth+net.area+c.dist+Gy, data = egg) 
summary(egg.Gy.lm) #0,9039 p-value presque 0

# On va transformer la réponse avec Box cox 
library(MASS)
egg.boxcox= boxcox( egg_start, lambda= seq( -2, 2, 0.05 ) )
egg2.boxcox= boxcox( egg_start, lambda= seq( -0.5, 1, 0.05 ) )
#Lambda est proche de 0 et inférieur à 0.5
#Lambda max : 
(lambda <-egg.boxcox$x[which.max(egg.boxcox$y)])
(lambda_precis <-egg2.boxcox$x[which.max(egg2.boxcox$y)])
#Lambda max est 0.1010101 pour egg.boxcox
#Lambda est 0.1060606 pour intervalle plus précise
#Test pour approximer avec log ou racine

eggLogReg = lm(log(egg.dens) ~ b.depth + lat + time + flow + s.depth +net.area + c.dist, data = egg)
summary(eggLogReg)

eggsqrtReg = lm(sqrt(egg.dens) ~ b.depth + lat + time + flow + s.depth +net.area + c.dist, data = egg)
summary(eggsqrtReg)

eggExpReg <- lm( I(egg.dens^lambda) ~ b.depth + lat + time + flow + s.depth +net.area + c.dist, data = egg)
summary(eggExpReg)

### On va comparer nos 3 possibilités avec QQplot et résidus
#plot(eggLogReg)
#plot(eggsqrtReg)
#plot(eggExpReg)
#racine rejet sans aucun doute

#Log et transfo pour lamba donnent résultats similaires
#On va continuer avec transformation log
#Examinons les distances cooks et potentiels 

potentiels.egg = hat( model.matrix( eggLogReg ) )
plot( potentiels.egg, pch=16, main="Potentiels: données maquereau" )
max(potentiels.egg) #0.3144777
which(potentiels.egg==max(potentiels.egg)) 
#312 a le plus haut potentiel

rstut0 <-rstudent(eggLogReg)
which(rstut0==max(rstut0)) #312
rstut0[312] #4.853776
T_critique <-qt(1-0.05/(2*368),361)
#312 outlier

#Distance de cook
cooks.egg = cooks.distance( eggLogReg )
plot( cooks.egg, ylab="Distances de Cook: Données maquereau", pch=16, cex=1)
library(olsrr)
plot(rstut0, ylab = "résidus studentisés")
abline(h=0)
ols_plot_cooksd_bar(eggLogReg)
cooks.egg[312] 
#1.27 Valeur élevée -> 312 est influente

#Résidus de student de la régression
resid.st <-rstudent(eggExpReg)
#Examinons la régression sans l'obs. 312
maquereau2 = egg[-312,]
egglm2 = lm(egg.dens ~ b.depth + lat + time + flow + s.depth + net.area + c.dist, data = maquereau2)
summary(egglm2)
#plot(egglm2)
resid.st2 <-rstudent(egglm2)
plot(resid.st2)
abline(h=0)

#####Boxcox 2
egg3.boxcox = boxcox( egglm2)
#encore très similaire
(lambda3 <-egg3.boxcox$x[which.max(egg3.boxcox$y)])
#lambda2 = 0.1010101
#Tranformation est encore valide

#Réexamen 
eggLogReg2 = lm(log(egg.dens) ~ b.depth + lat + time + flow + s.depth + net.area + c.dist, data = maquereau2)
summary(eggLogReg2)

eggsqrtReg2 = lm(sqrt(egg.dens) ~  b.depth + lat + time + flow + s.depth + net.area + c.dist, data = maquereau2)
summary(eggsqrtReg2)

eggExpReg2 <- lm( I(egg.dens^lambda3) ~  b.depth + lat + time + flow + s.depth + net.area + c.dist, data = maquereau2)
summary(eggExpReg2)


#plot(eggExpReg2)
#plot(eggLogReg2)
#plot(resid(eggExpReg2))
#abline(h=0)
#plot(resid(eggLogReg2))
#abline(h=0)
#on continue avec transformation y avec log
potentiels2 = hat( model.matrix( eggLogReg2 ) )
plot( potentiels2, pch=16, main="Potentiels: données maquereau sans obs. 312" )
cooks.egg2 = cooks.distance( eggLogReg2)
plot( cooks.egg2, ylab="Distances de Cook: Données maquereau sans obs.312", pch=16, cex=1)

cook.eggLog2 = cooks.distance( eggLogReg2)
max(cook.eggLog2) #0.066
ols_plot_cooksd_bar(eggLogReg2)
#Petites valeurs de distances de cook et potentiel de 0.15 
summary(eggLogReg2)

rstut1 <-rstudent(eggLogReg2)
which(abs(rstut1)==max(abs(rstut1))) #312
rstut0[95] #3.418675
T_critique <-qt(1-0.05/(2*368),360) #3.85 
###on  a pas de valeurs aberrantes pour linstant 

############################On va aller backward########################
summary(eggLogReg2)
#F inférieur à Fout ou regarder valeur- p avec tvalue
#On enlèeve c.dist car grande valeur-p et petite t value
#T-stat =  -1.494 et valeur- p  = 0.136100 (pas significative)
eggLog3 <-lm(formula = log(egg.dens) ~ b.depth + lat + time + flow + s.depth + net.area , data = maquereau2)
summary(eggLog3)

#On regarde s'il y a des données aberrantes 
T_critique2 <-qt(1-0.05/(2*368),361)
rstut2 <-rstudent(eggLog3)
which(abs(rstut2)==max(abs(rstut2)))
rstut2[95] #3.42
#pas de valeurs aberrantes
cooks.eggLog3 = cooks.distance(eggLog3)
plot( cooks.eggLog3, ylab="Distances de Cook: Données maquereau sans obs.312", pch=16, cex=1)
which(cooks.eggLog3==max(cooks.eggLog3))

#On enlèeve b.depth 

eggLog4 <-lm(formula = log(egg.dens) ~ lat + time + flow + s.depth + net.area , data = maquereau2)
summary(eggLog4)
#On garde ce modèle

#plot(eggLog4)
plot(resid(eggLog4))
abline(h=0)
#plot(eggLog4)
potentiels3 = hat( model.matrix( eggLog4 ) )
plot( potentiels3, pch=16, main="Potentiels: données maquereau 2" )
cooks.egg3 = cooks.distance(eggLog4)
plot( cooks.egg3, ylab="Distances de Cook: Données maquereau sans obs.312", pch=16, cex=1)
max(potentiels3) #0.109
which(potentiels3==max(potentiels3)) #346
T_critique3 <-qt(1-0.05/(2*368),362)
rstut3 <-rstudent(eggLog4)
which(abs(rstut3)==max(abs(rstut3)))
rstut3[95] #3.489



#Examionons la corrélation entre prédicteurs
library(car)
library(quantmod)
library(MASS)
library(corrplot)
vif(eggLog4)

#test pour variance dépend moyenne
ri = rstandard(eggLog4)
ei= resid(eggLog4)

SSE = ((1.182* 1.182)* 362)
sigma2 = SSE/nrow(maquereau2)
ui = ei*ei/sigma2
zi = fitted(eggLog4)
uy.lm = lm(ui  ~ zi)
summary(uy.lm)
SSreg= anova(uy.lm)$"Sum Sq"[1]

#test: valeur-p
#Ha suivrait la moyenne , on rjette H0
cat( "Valeurs prédites: ", SSreg/2, 1 - pchisq( SSreg/2,anova(uy.lm)$Df[1] ) )
#6.41254 0.01133172
plot(zi, ri^2, xlab="valeurs prédites", ylab="résidus", pch=16 )
abline(h=0)

library(leaps)
predictor_matrix= cbind(egg$b.depth,egg$lat,egg$time,egg$flow,egg$s.depth,egg$net.area,egg$c.dist)
Bestmodel = leaps(predictor_matrix, egg$egg.dens, int= TRUE, method = "Cp", 
                  names = c("b.depth","lat","time","flow","s.depth","net.area","c.dist"))
Bestmodel
plot(Bestmodel$size,Bestmodel$Cp, xlab = "taille du modèle", ylab= "Cp", pch=20)
abline(0,1)

# avec 6 variables , on est sur la droite et plus petite Cp (donc 5 prédicteurs + constante)


library(stats)
AIC(eggLog4)  #Best AIC 1175.656
AIC(egg_start) #4244.971
# variable toutes significatives 
#VIF pas de colinéarité



####On va analyser variance : l'hétéroscédatiscité ?

library(lmtest)
library(zoo)
bptest(eggLog4)
library(performance)
check_heteroscedasticity(eggLog4)
#Heteroscedasticity (non-constant error variance) detected (p = 0.011).
library(car)
ncvTest(eggLog4) #on rejette h0 

#moindres carrées pondérés

#à faire
#examiner les variables prédicteurs
#
# poids possibles : s.depth^4 ou net.area^3
anova(eggLog4)
#on voit que s.depth ou net.area a des mean sq plus élavuer, ils sont plus responsables de la variance élevée

#sachant que la réponse est logaithmique ou inversement prédicteurs suivent expo, on voit approximer par un
#poids exposant cubique

# Testons
eggLog_weights<-lm(formula = log(egg.dens) ~ lat + time + flow + s.depth + net.area , weights = 1/net.area^3, data = maquereau2)
plot(eggLog_weights)
summary(eggLog_weights)
#On regarde si variance constante
ncvTest(eggLog_weights)
check_heteroscedasticity(eggLog_weights) #aucun problème
vif(eggLog_weights) #aucun problème de colinéarité

#On va réaxaminer les résidus pour régression pondérée
plot(resid(eggLog_weights))
abline(h=0)
potentiels4 = hat( model.matrix( eggLog_weights) )
plot( potentiels4, pch=16, main="Potentiels: données maquereau 2" )
cooks.eggw4 = cooks.distance(eggLog_weights)
plot( cooks.eggw4, ylab="Distances de Cook: Données maquereau sans obs.312", pch=16, cex=1)
ols_plot_cooksd_bar(eggLog_weights)
max(potentiels4) #0.109
which(potentiels4==max(potentiels4)) #346
T_critiquew4 <-qt(1-0.05/(2*368),362) #3.86
rstut4 <-rstudent(eggLog_weights)
which(abs(rstut4)==max(abs(rstut4)))
which(abs(rstut4)>3) #95,118,131,136
which(abs(rstut4)>4) #95 seulement
rstut4[95] #4.613817 

#95 est outlier
cooks.eggw4[95]  # pas influent

#est-ce que 95 doit être enlever si aberrant pas influent?

#Est-ce qu'il y a d'autres aberrants?
#rstut4[118]
#rstut4[131]
#rstut[136]
#Non 

#on enlève 95 seulement

###### on essaie avec poids de s.depth
eggLog_weights_sdepthpoids<-lm(formula = log(egg.dens) ~ lat + time + flow + s.depth + net.area , weights = 1/s.depth^4, data = maquereau2)
summary(eggLog_weights_sdepthpoids)
#plot(eggLog_weights_sdepthpoid)  
## Problèeme résidus
#On regarde si variance constante
ncvTest(eggLog_weights_sdepthpoids)
check_heteroscedasticity(eggLog_weights_sdepthpoids)
cooks.egg_wEssaiSdpeth = cooks.distance(eggLog_weights_sdepthpoids) ## Distance de cook 25 
#solution non retenue


essai = egg[c(-95,-312),]

#on refait régression pondérée a
eggLog_weights2<-lm(formula = log(egg.dens) ~ lat + time + flow + s.depth + net.area , weights = 1/net.area^3, data = essai)
summary(eggLog_weights2)
plot(eggLog_weights2)
check_heteroscedasticity(eggLog_weights2) #nest pas homoscédastique, autre régression est mieux
#Vérification aucun résidus 
plot(resid(eggLog_weights2))
abline(h=0)
potentiels_w2= hat( model.matrix( eggLog_weights2) )
plot( potentiels_w2, pch=16, main="Potentiels: données maquereau 3" )
cooks.egg_w2 = cooks.distance(eggLog_weights2)
plot( cooks.egg_w2, ylab="Distances de Cook: Données maquereau sans obs.312 & 95", pch=16, cex=1)
max(potentiels_w2) #0.109
which(potentiels_w2==max(potentiels_w2)) #345
T_critique_poids2 <-qt(1-0.05/(2*367),361) #3.856
rstut_w2<-rstudent(eggLog_weights2)
which(abs(rstut_w2)==max(abs(rstut_w2)))
rstut_w2[135]
cooks.egg_w2[135]
#135 pas infleunnt
#inférieur èa  t-critique



#Comparaison AIC
AIC(eggLog4)  #Best AIC 1175.656
AIC(egg_start)
AIC(eggLog_weights)
AIC(eggLog_weights2)

#à décider compromis entre régression , variance, résidus et AIC
