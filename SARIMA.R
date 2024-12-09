source("./cod/Config.R")
source("./cod/ETL.R")

# Tests aumentados de Dickey-Fuller de todas las series
adf_testGA1 <- adf.test(GA1, alternative = "stationary")
adf_testGA1
# p-valor = menor a 0.01 -> Evidencia de estacionariedad

# Variables predictoras y objetivo
h <- 24
N <- length(GA1)
Y24 <- GA1[(N - 23):N]

# Ciclo para realizar nowcasting bajo una ventana de los ultimos 24 meses
nowcasts <- numeric(h)
for(hi  in 1:h){
    # Evaluacion de modelo SARIMA
    sarima_model <- Arima(GA1[1:(N + hi - h - 1)], order = c(0, 0, 0), seasonal = list(order = c(0, 1, 1), period = 12))
    summary(sarima_model)# Revisamos los residuales del modelo SARIMA ajustado

    # Realizamos nowcasting con SARIMA
    GA1_presente_SARIMA <- forecast(sarima_model, h = 1)
    nowcasts[hi] <- GA1_presente_SARIMA$mean
    print(hi)
}

# Calculamos accuracy de los valores predichos por sarima 
accuracy_SARIMA <- accuracy(nowcasts, Y24)
accuracy_SARIMA

png("./images/Resultados/PredReal_SARIMA.png")
plot(as.vector(nowcasts), Y24, type = "p", pch = 19,
main = "Comparación GA1 con modelo SARIMA",
xlab = "GA1 predicho",
ylab = "GA1 real")
abline(0, 1, col = "red")
dev.off()

write.csv(data.frame(nowSARIMA = as.vector(nowcasts), Obs = Y24), file = "./data/nowSARIMA.csv", row.names = FALSE)
