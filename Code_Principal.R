# chargement des paquets
charger_pkgs <- function(...){
  pkgs <- as.character(match.call())[-1]
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
  if(length(new.pkgs)) install.packages(new.pkgs)
  error <- try(expr = {lapply(pkgs,function(x){library(x,character.only=TRUE)})},silent = FALSE)
}
charger_pkgs("magrittr","ESG")
rm(charger_pkgs)


# récupération des zero-coupons : 
data(ZC)

# Fonction de gse : 
build_gse <- function(ZC = ZC,
                      base.horizon = 5, # horizon de projection 
                      base.nScenarios = 1000, # nombre de scénarios
                      rt.vol = .1, # Volatilité du taux courts dans le modèle de vasicek
                      rt.k = 2, # paramètre k du vasicek
                      s.vol = .1, # Volatilité du taux courts dans le modèle de vasicek pour les actions
                      s.k = 2, # paramètre k du vasicek actions
                      s.volStock = .2, # Volatilité des actions (sigma de blakc-sholes)
                      s.stock0 = 100, # Valeur initiale de l'action
                      s.rho=.5 # correlation entre l'action et les taux d'intérets.
                      ){
  # Création des scénarios : 
  new("Scenarios") %>% 
    setParamsBaseScenarios(     # Paramètres de base
      horizon = base.horizon,              # horizon de projection 
      nScenarios = base.nScenarios         # nombre de scénarios
    ) %>% 
    setRiskParamsScenariosrt(   # Taux courts : Modèle HJM avec un paramétrage du sigma Hull-White (= vaiscek généralisé)
      vol = rt.vol,                 # Volatilité du taux courts dans le modèle de vasicek
      k = rt.k                     # paramètre k du vasicek
    ) %>% 
    setRiskParamsScenariosS(    # Actions : black-sholes avec taux d'intéret stochastique (vaiscek) (correpond au taux court)
      vol = s.vol,                 # Volatilité du taux courts dans le modèle de vasicek
      k = s.k,                    # paramètre k du vasicek
      volStock = s.volStock,            # Volatilité des actions (sigma de blakc-sholes)
      stock0 = s.stock0,             # Valeur initiale de l'action
      rho = s.rho                   # correlation entre l'action et les taux d'intérets.
    ) %>% 
    setForwardRates(            # Récupération de la courbe des taux
      ZC,
      horizon=base.horizon
    ) %>% 
    setZCRates(                 # Inclusion des taux courts
      ZC,
      horizon=base.horizon
    ) %>%                       
    customPathsGeneration(      # Génération des trajectoires de taux court
      type="shortRate"
    ) %>% 
    customPathsGeneration(      # Génération des trajectoire action.
      type="stock"
    ) %>% 
    return
}

objScenario <- build_gse(ZC = ZC)
# Test de martingalité :
MartingaleTest(objScenario)



#--------------------------- Fonction de calcul des flux -----------------
calculFlux <- function(objScenario,
                       txStructurel,
                       txConjoncturel)
{
  #Projection des flux
  PM <- matrix(data=1,nrow=nrow(objScenario@shortRatePaths),ncol=ncol(objScenario@shortRatePaths))
  flux <- matrix(data=0,nrow=nrow(objScenario@shortRatePaths),ncol=ncol(objScenario@shortRatePaths))
  actu <- matrix(data=1,nrow=nrow(objScenario@shortRatePaths),ncol=ncol(objScenario@shortRatePaths))
  T <- ncol(objScenario@shortRatePaths)-1
  
  for (t in 1:T)
  {
    if (t==T)
    {
      tauxRachat <- 1
    }else
    {
      tauxRachat <- txStructurel+txConjoncturel*(objScenario@stockPaths[,t]<1)
    }
    PM[,t+1] <- PM[,t]*objScenario@stockPaths[,t+1]/objScenario@stockPaths[,t]
    flux[,t+1] <- PM[,t+1]*tauxRachat
    PM[,t+1] <- PM[,t+1]-flux[,t+1]
    actu[,t+1] <- actu[,t]*exp(-objScenario@shortRatePaths[,t+1])
  }
  return(list(PM=PM,flux=flux,actu=actu))
}

calculFlux(objScenario,0.3,0.6)

calculFlux(objScenario,0.3,0.6) %>% {.$flux * .$actu} %>%
  rowSums %>%
  mean %>% {print(paste0("Le Be est égal à ",.))}


