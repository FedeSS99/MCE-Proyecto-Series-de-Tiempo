library(tseries)
library(forecast)

source("./Functions_DFM.R")

library(ggplot2)
library(ggfortify)

datosGA1 <- read.csv("./data/GA1.csv")
datosSIAP <- read.csv("./data/IVF_SIAP.csv")

# intervalo de tiempo que comparten ene2015-dic2023
GA1 <- ts(datosGA1[265:372, 2], start = c(2015, 1), frequency = 12)
IVF <- ts(datosSIAP[1:108, 2], start = c(2015, 1), frequency = 12)

# Revisamos las correlaciones entre la var objetivo y las predictoras
cor(GA1, IVF)

# Obtenemos estimaciones de los factores 
rhat <- 1
demean <- 2
