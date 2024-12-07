library(tseries)
library(forecast)

library(ggplot2)
library(ggfortify)

datosGA1 <- read.csv("./data/GA1.csv")
datosSIAP <- read.csv("./data/IVF_SIAP.csv")

# intervalo de tiempo que comparten ene2015-dic2023
GA1 <- ts(datosGA1[265:372, 2], start = 2015, frequency = 12)
IVF <- ts(datosSIAP[1:108, 2], start = 2015, frequency = 12)

# Visualizacion de GA1 y SIAP entre ene2015-dic2023
autoplot(GA1, ts.colour = "blue") + ggtitle("GA1")
autoplot(IVF, ts.colour = "black") + ggtitle("SIAP")

# Graficas ACF Y PCF de series originales
par(mfrow = c(2, 1))
acf(GA1)
pacf(GA1)

par(mfrow = c(2, 1))
acf(IVF)
pacf(IVF)

# Tests aumentados de Dickey-Fuller en ambas series
adf_testGA1comp <- adf.test(GA1, alternative = "stationary")
adf_testGA1comp
# p-valor = menor a 0.01 -> Evidencia de estacionariedad

adf_testSIAPcomp <- adf.test(IVF, alternative = "stationary")
adf_testSIAPcomp
# p-valor = 0.0153 -> Evidencia de estacionariedad

# Variables predictoras y objetivo
N <- length(GA1)
X <- IVF[1:(N - 1)]
Y <- GA1[1:(N - 1)]
Xpresente <- IVF[N]
Ypresente <- GA1[N]

# Ajuste de modelo ARIMA
arima_model <- auto.arima(Y, xreg = X, max.p = 12, max.q = 12, stepwise = FALSE, nmodels = 256, approximation = FALSE, stationary = TRUE, seasonal = FALSE, parallel = TRUE, num.cores = NULL)
summary(arima_model)
# Revisamos los residuales del modelo ARIMA ajustado
checkresiduals(arima_model)

# Ajuste de modelo SARIMA
sarima_model <- auto.arima(Y, xreg = X, max.p = 12, max.q = 12, max.P = 12, max.Q = 12, stepwise = FALSE, nmodels = 256, approximation = FALSE, stationary = TRUE, seasonal = TRUE, parallel = TRUE, num.cores = NULL)
summary(sarima_model)
# Revisamos los residuales del modelo SARIMA ajustado
checkresiduals(sarima_model)

# Realizamos nowcasting con ARIMA
GA1_presente_ARIMA <- forecast(arima_model, xreg = Xpresente, h = 1)
GA1_presente_ARIMA
par(mfrow = c(1,1))
plot(GA1_presente_ARIMA)

# Realizamos nowcasting con SARIMA
GA1_presente_SARIMA <- forecast(sarima_model, xreg = Xpresente, h = 1)
GA1_presente_SARIMA
par(mfrow = c(1,1))
plot(GA1_presente_SARIMA)

# Calculamos accuracy del valor predicho por SARIMA con el real
accuracy_SARIMA <- accuracy(GA1_presente_SARIMA, Ypresente)
accuracy_SARIMA
