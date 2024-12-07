library(tseries)
library(forecast)

library(ggplot2)
library(ggfortify)

datosGA1 <- read.csv("./data/GA1.csv")
datosSIAP <- read.csv("./data/IVF_SIAP.csv")

# intervalo de tiempo que comparten ene2015-dic2023
GA1_comp <- ts(datosGA1[265:372, ], start = 2015, frequency = 12)
SIAP_comp <- ts(datosSIAP[1:108, ], start = 2015, frequency = 12)

# Visualizacion de GA1 y SIAP entre ene2015-dic2023
autoplot(GA1_comp[, 2], ts.colour = "blue") + ggtitle("GA1")
autoplot(SIAP_comp[, 2], ts.colour = "black") + ggtitle("SIAP")

# Buscamos un valor d optimo
for(d in 1:12){
    DiffGA1 <- diff(GA1_comp[,2], differences = d)
    acf(DiffGA1, lag.max = 60)
}
# Observamos que no hay cambios significativos entre los diferentes ordenes de diferenciacion
# por lo que d = 0

SelecDiff_GA1 <- diff(GA1_comp[, 2], differences = 1)

# Observamos la PACF para la serie diferenciada
pacf(SelecDiff_GA1, lag.max = 60)

# Observamos la PACF para la serie diferenciada
acf(SelecDiff_GA1, lag.max = 60)
# Vemos que no hay un corte subito y en camibo hay una caída exponencial
# por lo que asumimos que q = 0

arima_model <- Arima(GA1_comp[, 2], xreg = SIAP_comp[,2], order = c(6, 1, 0), seasonal = c(0,0,0))
summary(arima_model)