# TD1 Générateur de scénarios économiques####
#############################################

# Pas d'évaluation dans ce cours !
install.packages('ESG')
library("ESG")

#Objectif: Générer des scénarios économiques
data("ZC")

A <- rStock(horizon=10, nScenarios=1000, ZC=ZC, vol=.2, k=2, volStock=.2, stock0=100, rho=.5)
summary(A)

#Changer horizons et nombre de simulation
matplot(t(A$stockPaths),type='l')

moy=c(0)
par(mfrow=c(2,2))
i=1
for(i in 1:4){
  
#hist(t(A$stockPaths))
matplot(t(A$stockPaths),type='l')
moy=c(moy,mean(A$stockPaths))
}

moy

#Changer horizons et nombre de simulation
rt <- rShortRate(horizon=15, nScenarios=10, ZC=ZC, vol=.1, k=2)
matplot(t(rt),type='l')


#Simuler tous les indices en meme temps
I<-rAllRisksFactors(horizon=5, nScenarios=10, ZC, vol=.1, k=2, volStock=.2, stock0=100, rho=.5, volRealEstate=.15,
                 realEstate0=50, eta=.05,liquiditySpread0=.01, defaultSpread0=.01,volDefault=.2,alpha=.1,beta=1)

par(mfrow=c(2,2))

matplot(t(I$s), type='l')
matplot(t(I$realEstate), type='l')
matplot(t(I$liquiditySpread), type='l')
matplot(t(I$defaultSpread), type='l')

#2ème partie
k=0.12
sTaux=0.05
sUC=0.16
H=40
nSimulations=10
tauxRachatS=0.03
tauxRachatC=0.06


#2ème partie

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
  res <- list()
  res[["PM"]] <- PM
  res[["flux"]] <- flux
  res[["actu"]] <- actu
  calculFlux <- (res)
}


k=0.12
sTaux=0.05
sUC=0.16
H=40
nSimulations=10000
tauxRachatS=0.03
tauxRachatC=0.06

trajTaux <- rStock(horizon=H, nScenario=nSimulations, ZC=ZC, vol=sTaux, k=k,volStock=sUC, stock0=1, rho=0.5)
trajUC <- rStock(horizon=H, nScenario=nSimulations, ZC=ZC, vol=sUC, k=k,volStock=sUC, stock0=1, rho=0.5)
summary(trajTaux)
calculFlux <- calculFlux(trajTaux$shortRatePaths, trajTaux$stockPaths ,tauxRachatS,tauxRachatC)
	
calculFluxFutur <- calculFlux$flux*calculFlux$actu

BEmpirique <- sum(calculFluxFutur )/nSimulations
BEmpirique 


