library(tseries)
library(forecast)

library(ggplot2)
library(ggfortify)

datosGA1 <- read.csv("./data/GA1.csv")
datosSIAP <- read.csv("./data/IVF_SIAP.csv")

# intervalo de tiempo que comparten ene2015-dic2023
GA1_comp <- ts(datosGA1[265:372,])
SIAP_comp <- ts(datosSIAP[1:108,])

# Visualizacion de GA1 y SIAP entre ene2015-dic2023
autoplot(GA1_comp[,2], ts.colour = "blue") + ggtitle("GA1")
autoplot(SIAP_comp[,2], ts.colour = "black") + ggtitle("SIAP")

# Graficas ACF Y PCF de series originales
par(mfrow = c(2, 1))
acf(GA1_comp[, 2])
pacf(GA1_comp[, 2])

par(mfrow = c(2, 1))
acf(SIAP_comp[, 2])
pacf(SIAP_comp[, 2])

# Transformacion logaritmica base 10 y Diferenciacion con lag = 12 (Posible overfitting)
tranGA1 <- diff(log10(GA1_comp[,2]), lag = 12)
tranSIAP <- diff(log10(SIAP_comp[,2]), lag = 12)

# Visualizamos las series bajo la transformacion y diferenciacion
autoplot(tranGA1, ts.colour = "blue") + ggtitle("GA1")
autoplot(tranSIAP, ts.colour = "black") + ggtitle("SIAP")

# Graficas ACF Y PCF de series transformadas
par(mfrow = c(2, 1))
acf(tranGA1)
pacf(tranSIAP)

par(mfrow = c(2, 1))
acf(tranGA1)
pacf(tranSIAP)

# Tests aumentados de Dickey-Fuller en ambas series
adf_testAG1 <- adf.test(tranGA1, alternative = "stationary")
adf_testAG1
# p-valor = 0.02 < 0.05 -> Evidencia de estacionariedad

adf_testSIAP <- adf.test(tranSIAP, alternative = "stationary")
adf_testSIAP
# p-valor = 0.045 < 0.05 -> Evidencia de estacionariedad (dudoso)

# Variables predictoras y objetivo
N <- nrow(GA1_comp)
X <- SIAP_comp[1:(N - 1), 2]
Y <- GA1_comp[1:(N - 1), 2]
Xpresente <- SIAP_comp[N, 2]
Ypresente <- GA1_comp[N, 2]

# Ajuste de modelo SARIMA
sarima_model <- Arima(Y, order = c(12, 0, 0), xreg = X)
summary(sarima_model)
# Revisamos los residuales del modelo SARIMA ajustado
checkresiduals(sarima_model)

GA1_presente_SARIMA <- forecast(sarima_model, xreg = Xpresente, h = 1)
GA1_presente_SARIMA
par(mfrow = c(1,1))
plot(GA1_presente_SARIMA)

# Calculamos accuracy del valor predicho por SARIMA con el real
accuracy_SARIMA <- accuracy(GA1_presente_SARIMA, Ypresente)
accuracy_SARIMA
