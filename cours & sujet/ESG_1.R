library("ESG")

data("ZC")

objScenario <- new("Scenarios")

objScenario <- setParamsBaseScenarios(objScenario
                                      , horizon = 1
                                      , nScenarios = 100)




#Sans precision sur le parametrage, simuler et observer des valeurs action
simulStock <- rStock(horizon = 5
                     , nScenarios = 10
                     , ZC = ZC
                     , vol = 0.1
                     , k = 3
                     , volStock = 0.1
                     , stock0 = 100
                     , rho = -0.3
                     )

matplot(t(simulStock$stockPaths), type = "l")

simulShort <- rShortRate(horizon = 5
                     , nScenarios = 10
                     , ZC = ZC
                     , vol = 0.1
                     , k = 3)

matplot(t(simulShort), type = "l")

simulEstate <- rRealEstate(horizon = 5
                         , nScenarios = 10
                         , ZC = ZC
                         , vol = 0.1
                         , k = 3
                         , volRealEstate = 0.2
                         , realEstate0 = 10)

matplot(t(simulEstate$realEstatePaths), type = "l")

allIndexes <- ESG::rAllRisksFactors(
  horizon = 5
  , nScenarios = 10
  , ZC = ZC
  , vol = 0.1
  , k = 3
  , volRealEstate = 0.2
  , realEstate0 = 10
  , volStock = 0.1
  , stock0 = 100
  , rho = -0.3
  , eta = 0.4
  , volDefault = 0.3
  , defaultSpread0 = 1
  , liquiditySpread0 = 1
  , alpha = 1
  , beta = 1
  )

timeSeq <- 1:30
sapply(timeSeq, 
       function(x)
       {
         time1 <- Sys.time()
         simu <- rShortRate(horizon = x, nScenarios = 50000, ZC = ZC, vol = 0.3, k = 3)
         time2 <- Sys.time()
         return(as.numeric(difftime(time2, time1)))
       }
       ) %>% plot(type = "l")

nSims <- seq(10, 100000, length.out = 200)
BEConv <- sapply(nSims, )
nSim = 100000
traj <- rStock(horizon = 40
               , nScenarios = nSim
               , ZC = ZC
               , vol = 0.05
               , k = 0.12
               , volStock = 0.16
               , stock0 = 1
               , rho = 0.5)

flux <- calculFlux(scenariosTaux = traj$shortRatePaths
                   , scenariosUC = traj$stockPaths
                   , txStructurel = 0.03 
                   , txConjoncturel = 0.06)

#Flux qui vont être lié a la loi comportementale de nos projections; tout ça c'est des flux qui vont être perdus
flux$flux %>% t %>% plot
flux$PM %>% t %>% plot
flux$actu %>% t %>% plot

flux$flux * flux$actu -> fluxActu

BEmpirique <- sum(fluxActu)/nSim
BEmpirique
