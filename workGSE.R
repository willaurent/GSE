# chargement des paquets
charger_pkgs <- function(...){
  pkgs <- as.character(match.call())[-1]
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new.pkgs)) install.packages(new.pkgs)
  error <- try(expr = {lapply(pkgs,function(x){library(x,character.only=TRUE)})},silent = FALSE)
}
charger_pkgs(
  "magrittr",
  "tidyverse",
  "purrr",
  "ESG"
)
rm(charger_pkgs)

# probablement a regarder dans la doc : 
# Tuto du packages.
# rAllRisksFactors
# rAssetDistribution
# rDefaultSpear
# rLiquiditySpread
# rRealEstate
# rShortRate
# rStock
# MartingaleTest()
# ParamScenarios class


# Import d'une courbe de prix zéro-coupon
data("ZC")
rStock(horizon = 15, # Horison de projection.
       nScenarios = 40, # Nombre de scénarii
       ZC = ZC, # les taux zéro-coupon en inputs
       vol=0.1, # valatilité du taux court
       k=2, # k for rates in vasicek models
       volStock=.2, # volatilité de l'action
       stock0 = 100, # valeur initiale de l'actoin
       rho=.5) %>% # Correlation entre taux court et action
       {matplot(t(.$stockPaths),type="l",xlab="Horizon",ylab="Action")}



# Sans précision sur le paramétrage, simuler et observe ddu taux courts
rShortRate(horizon = 15,
           nScenarios = 40,
           ZC = ZC,
           vol = 0.2,
           k = 10) %>%
           {matplot(t(.),type="l",xlab="Horizon",ylab="taux court")}


# Sans précision sur le paramétrage, simuler et observe des valerus immobilier
rRealEstate(horizon = 15,
            nScenarios = 40,
            ZC = ZC,
            vol = 0.2,
            k = 10,
            volRealEstate = 0.2,
            realEstate0 = 100) %>%
            {matplot(t(.$realEstatePaths),type="l",xlab="Horizon",ylab="Immobilier")}






# simuler et observe tous les indices en simulténa"
par(mfrow=c(3,2))
rAllRisksFactors(horizon = 10,
                 nScenarios = 100,
                 ZC = ZC,
                 vol = 0.2,
                 k = 10,
                 volStock = 0.1,
                 stock0 = 100,
                 rho = 0.5,
                 volRealEstate = 0.1,
                 realEstate0 = 100,
                 eta = 0.4,
                 liquiditySpread0 = 100,
                 defaultSpread0 = 200,
                 volDefault = 1.6,
                 alpha = 1,
                 beta = 2) %>%
                 {map(1:length(.), function(x) matplot(t(.[[x]]),type="l",xlab="horison",ylab=names(.)[x],main=names(.)[x]))}









# Msreur le temps d'éxécution de la simulation : 

Temps <- function(horizon,nScenarios,type.of.time=1){
  time1<-proc.time() 
  rAllRisksFactors(horizon = horizon,
                   nScenarios = nScenarios,
                   ZC = ZC,
                   vol = 0.2,
                   k = 1,
                   volStock = 0.1,
                   stock0 = 100,
                   rho = 0.5,
                   volRealEstate = 0.1,
                   realEstate0 = 100,
                   eta = 0.4,
                   liquiditySpread0 = 100,
                   defaultSpread0 = 200,
                   volDefault = 1.6,
                   alpha = 1,
                   beta = 2)
  
  return((proc.time()- time1)[type.of.time])
}



calc.temps <- function(horizon=c(1,5),scenario=c(1,10,100,1000)){
  
  kikou <- matrix(NA,nrow=length(scenario),ncol=length(horizon))
  
  for(i in 1:length(horizon)){
    for(j in 1:length(scenario)){
      kikou[j,i] <- Temps(horizon[i],scenario[j])
    }
  }
  
  return(list(x=horizon,y=scenario,temps=t(kikou)))
}


calc.temps(horizon = seq(1,50,by=5),
           scenario = 10^(1:4)) %$% 
  persp(x,y,temps)






### Calcul de BE sur R : 

# Import de la conftion CalculFlux 

# Elle permet de calculer les flux futurs en tenant compte : 
#     D'un faceur de revalorisation obtenu a l'aide des scénarios UC
#   D'un facteur d'actu obtenu a l'aide des scénarios
#   de paramètres comportementaux ( rachat structurel et conjoncturel) 
# 
# 
# 


# parmaètres comportementaux : taux a fixer : 3%

# facteur d'actualisation : 




#--------------------------- Fonction de calcul des flux -----------------
calculFlux <- function(scenariosTaux,scenariosUC,txStructurel,txConjoncturel)
{
  #Projection des flux
  PM <- matrix(data=1,nrow=nrow(scenariosTaux),ncol=ncol(scenariosTaux))
  flux <- matrix(data=0,nrow=nrow(scenariosTaux),ncol=ncol(scenariosTaux))
  actu <- matrix(data=1,nrow=nrow(scenariosTaux),ncol=ncol(scenariosTaux))
  T <- ncol(scenariosTaux)-1
  
  for (t in 1:T)
  {
    if (t==T)
    {
      tauxRachat <- 1
    }else
    {
      tauxRachat <- txStructurel+txConjoncturel*(scenariosUC[,t]<1)
    }
    PM[,t+1] <- PM[,t]*scenariosUC[,t+1]/scenariosUC[,t]
    flux[,t+1] <- PM[,t+1]*tauxRachat
    PM[,t+1] <- PM[,t+1]-flux[,t+1]
    actu[,t+1] <- actu[,t]*exp(-scenariosTaux[,t+1])
  }
  return(list(PM=PM,flux=flux,actu=actu))
}


par(mfrow=c(1,3))
rStock(horizon = 15, # Horison de projection.
       nScenarios = 40, # Nombre de scénarii
       ZC = ZC, # les taux zéro-coupon en inputs
       vol=0.1, # valatilité du taux court
       k=2, # k for rates in vasicek models
       volStock=.2, # volatilité de l'action
       stock0 = 100, # valeur initiale de l'actoin
       rho=.5) %$% 
  calculFlux(shortRatePaths,stockPaths,0.3,0.6)  %T>%
  {map(1:length(.), 
       function(x) matplot(t(.[[x]]),
                           type="l",
                           xlab="horison",
                           ylab=names(.)[x],
                           main=names(.)[x]))
  } %>% 
  {.$flux * .$actu} %>%
  rowSums %>%
  mean %>% {print(paste0("Le Be est égal à ",.))}


Be <- c(seq(1,500,by=1),seq(from=500,to=1000,by=2),seq(from=1000,to=2000,by=10),seq(2000,10000,by=100))
x= Be
for(i in 1:length(Be)){
  rStock(horizon = 15, # Horison de projection.
         nScenarios = Be[i], # Nombre de scénarii
         ZC = ZC, # les taux zéro-coupon en inputs
         vol=0.1, # valatilité du taux court
         k=2, # k for rates in vasicek models
         volStock=.2, # volatilité de l'action
         stock0 = 100, # valeur initiale de l'actoin
         rho=.5) %$% 
    calculFlux(shortRatePaths,stockPaths,0.3,0.6)  %>% 
    {.$flux * .$actu} %>%
    rowSums %>%
    quantile(probs=0.995) -> 
    Be[i] 
}

# graphe de convergence du BE : 
plot(x,Be,xlab="Nombre de simulations",ylab="Best estimate empirique")
data.frame(x=x,y=Be) %>% 
  ggplot(aes(x,y)) +
  geom_point()+
  geom_smooth()+
  ggtitle("Be empirique en fonction du nombre de simulation")+
  labs(x = "Nombre de simu", y ="Be empirique")





# prendre plusieurs Be avc des paam différents..; 

# Autre application projeter sur la durée du lan stratégique orsa les prix d'instruments financiers 
# mise ne place d'ajustement financiers dans un cadre ORSA

