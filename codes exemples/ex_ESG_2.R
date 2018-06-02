library(magrittr)
library(ESG)
library(nleqslv)

data("ZC")

inputData <- data.frame(
  maturity = c(1,1,1,exp(1), 1, 1, 1, 1, 1),
  strike = c(0.80, 0.90, 0.95, 0.97, 1, 1.025, 1.05, 1.1, 1.2),
  volatility = c(0.2519, 0.2244, 0.2137, 0.2090, 0.2047, 0.2006, 0.1967, 0.1891, 0.1756),
  rate = rep(0.0016, 9),
  dividende = rep(0.3329, 9)
)
S_0 <- 3146.43



inputData %>% print

# Exercice
# Calculer le prix d'un call en utilisant les données du marché
# Faire de même en utilisant un futur paramètre variable sigma
# calculer la somme des erreurs quadratiques (absolue et relative)
# Minimiser cette erreur en utilisant le solveur Excel
# Discrétisaer et créer 1000 scénario à partir du Black & Scholes
# Tester la mtg sur l'accroissement des rendements actions simulés

# Valeur d'un call de maturité tau = T - t est :
# Ct = exp(- q_(dividende) * (T-t) ) St * N(d1) - exp ( - r * (T-t) ) K * N(d2)
# d1 = ( Ln(S_t / K) + ((r - q + 0.5 sigma ²)(T-t))) / (sigma sqrt(T-t))
# d2 = d1 - sigma sqrt(T-t)


callPricing <- function(t,
                        maturity,
                        S_t,
                        q,
                        sigma,
                        K,
                        r)
{
  d1 = ( log( S_t / K) + ( (r - q + 0.5 * sigma ^ 2)*(maturity-t))  ) / (sigma*sqrt(maturity-t))
  d2 = d1 - sigma * sqrt(maturity-t)
  return(exp(-q * (maturity-t)) * S_t * pnorm(d1) - exp( - r * (maturity - t)) * K * pnorm(d2))
}
  
inputData$callPrice <- callPricing(t = rep(0, 9),
                                   maturity = inputData$maturity,
                                   S_t = rep(S_0, 9),
                                   q = inputData$dividende,
                                   sigma = inputData$volatility,
                                   K = inputData$strike*S_0,
                                   r = inputData$rate
)

inputData$callPrice %>% plot(y =  inputData$strike, type = "l")

calibrationFunction <- function(volatility, .inputData = inputData, .S_0 = S_0)
{
  marketPrice <- .inputData$callPrice
  
  # marketVolaility <- .inputData$volatility
  
  modelPrice <- callPricing(t = rep(0, 9),
                            maturity = .inputData$maturity,
                            S_t = rep(.S_0, 9),
                            q = .inputData$dividende,
                            sigma = volatility,
                            K = .inputData$strike*.S_0,
                            r = .inputData$rate)
  
  # impliedVolatilityModel <- NR(modelPrice)
  
  return (sum((marketPrice - modelPrice)^2))
}

# blackVolatility <- rootSolve::multiroot(f = calibrationFunction,0.02,maxiter = 100)

optim(method = "BFGS", par = 0.02, fn = calibrationFunction) -> optResults


inputData$callPriceBacktest <- callPricing(t = rep(0, 9),
                                   maturity = inputData$maturity,
                                   S_t = rep(S_0, 9),
                                   q = inputData$dividende,
                                   sigma = rep(optResults$par, 9),
                                   K = inputData$strike*S_0,
                                   r = inputData$rate
)

inputData$callPriceBacktest

inputData$relativeError <- (inputData$callPrice - inputData$callPriceBacktest)/inputData$callPrice

inputData$relativeError %>% plot




