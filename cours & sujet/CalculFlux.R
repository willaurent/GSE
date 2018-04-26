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
  return(res)
}