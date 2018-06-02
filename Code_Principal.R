################################################################ MISE EN PLACE DE LENVIRONNEMNET ####
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
ZC <- read.csv2(file="courbe_des_taux.csv")$EUR

################################################################ FONCTIONS UTILITAIRES ####
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


BE <- function(ZC = ZC,
               base.horizon = 5,
               base.nScenarios = 10000,
               rt.vol = c(.1,.2),
               rt.k = 2,
               s.vol = .1,
               s.k = 2,
               s.volStock = .2,
               s.stock0 = 100,
               s.rho=.5,
               txStructurel=0.3,
               txConjoncturel=0.1){
  df <- expand.grid(base.horizon=base.horizon,
                    base.nScenarios=base.nScenarios,
                    rt.vol=rt.vol,
                    rt.k=rt.k,
                    s.vol=s.vol,
                    s.k=s.k,
                    s.volStock=s.volStock,
                    s.stock0=s.stock0,
                    s.rho=s.rho,
                    txStructurel=txStructurel,
                    txConjoncturel=txConjoncturel)
  df$BE <- rep(NULL,nrow(df))
  for(i in 1:nrow(df)){
    build_gse(
      ZC = ZC,
      base.horizon = df$base.horizon[i],
      base.nScenarios = df$base.nScenarios[i],
      rt.vol = df$rt.vol[i],
      rt.k = df$rt.k[i] ,
      s.vol = df$s.vol[i],
      s.k = df$s.k[i],
      s.volStock = df$s.volStock[i],
      s.stock0 = df$s.stock0[i],
      s.rho=df$s.rho[i]
    ) %>%
      calculFlux(
        df$txStructurel[i],
        df$txConjoncturel[i]
      ) %>% 
      {.$flux * .$actu} %>%
      rowSums %>%
      mean -> df$BE[i]
  }
  return(df)
}


################################################################ TEMPS DE CONVERGENCE ####
objScenario <- build_gse(ZC = ZC,base.nScenarios = 10000)
MartingaleTest(objScenario)
rez <- calculFlux(objScenario,0.3,0.6)
BE <- rez %>% {.$flux * .$actu} %>% rowSums
BE %>% cummean %>% plot(type="l") %>% abline(h=1,col=2)
BE %>% mean
seuil=0.005
rev(which((BE %>% cummean <1+seuil)*(BE %>% cummean >1-seuil)==0))[1]




################################################################ MARTINGALITE ####
MartingaleTest(objScenario)
objScenario <- build_gse(ZC = ZC,
                         base.horizon = 5,
                         base.nScenarios = 10000,
                         rt.vol = .1,
                         rt.k = 2,
                         s.vol = .1,
                         s.k = 2,
                         s.volStock = .2,
                         s.stock0 = 100,
                         s.rho=.5)

################################################################ GRAPHIQUES DE SENSIBILITE ####
par(mfrow=c(3,3))
# BE(ZC=ZC,
#    base.horizon = seq(5,30,by=1)) %>% 
#    {plot(.$base.horizon,
#          .$BE,
#          main="BE en fonction de l'horizon du contrat"
#    )
#      abline(h=1)}
BE(ZC=ZC,
   rt.vol = seq(0,1,length.out=100)) %>% 
   {plot(.$rt.vol,
         .$BE,
         main="BE en fonction de la volatilité du taux court"
   )
     abline(h=1)}
BE(ZC=ZC,
   rt.k = seq(1,3,length.out=50)) %>% 
   {plot(.$rt.k,
         .$BE,
         main="BE en fonction du k du taux court"
   )
     abline(h=1)}
BE(ZC=ZC,
   s.vol = seq(0,1,length.out=100)) %>% 
   {plot(.$s.vol,
         .$BE,
         main="BE en fonction de la volatilité du taux court des actions"
   )
     abline(h=1)}
BE(ZC=ZC,
   s.k = seq(1,3,length.out=50)) %>% 
   {plot(.$s.k,
         .$BE,
         main="BE en fonction du k du taux court des actions"
   )
     abline(h=1)}
BE(ZC=ZC,
   s.volStock = seq(0,1,length.out=100)) %>% 
   {plot(.$s.volStock,
         .$BE,
         main="BE en fonction de la volatilité des actions"
   )
     abline(h=1)}
BE(ZC=ZC,
   s.stock0 = seq(10,200,length.out=100)) %>% 
   {plot(.$s.stock0,
         .$BE,
         main="BE en fonction la valeur initiale des actions"
   )
     abline(h=1)}
BE(ZC=ZC,
   s.rho=seq(0,1,length.out=100)) %>% 
   {plot(.$s.rho,
         .$BE,
         main="BE en fonction du rho des actions"
   )
     abline(h=1)}
BE(ZC=ZC,
   txConjoncturel=seq(0,1,length.out=100)) %>% 
   {plot(.$txConjoncturel,
         .$BE,
         main="BE en fonction du taux de rachat conjoncturel"
   )
     abline(h=1)}


BE(ZC=ZC,
   txStructurel=seq(0,1,length.out=100)) %>% 
   {plot(.$txStructurel,
         .$BE,
         main="BE en fonction du taux de rachat structurel"
   )
     abline(h=1)}



# BE(ZC=ZC,
#    base.horizon = 5,
#    base.nScenarios = c(10000),
#    rt.vol = .5,
#    rt.k = 2,
#    s.vol = seq(0,1,length.out=100),
#    s.k = 2,
#    s.volStock = .2,
#    s.stock0 = 100,
#    s.rho=.5,
#    txStructurel=seq(0,1,length.out=100),
#    txConjoncturel=0.1) %>% 
#    {plot(.$s.vol,
#          .$BE,
#          main="BE en fonction de la volatilité des actions"
#    )
#      abline(h=1)}



################################################################ The END ####