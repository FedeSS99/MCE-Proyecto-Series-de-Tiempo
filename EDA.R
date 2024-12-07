library(tseries)
library(ggplot2)
library(ggfortify)

datosGA1 <- read.csv("./data/GA1.csv")
datosIVF <- read.csv("./data/IVF_SIAP.csv")

# intervalo de tiempo que comparten ene2015-dic2023
GA1 <- ts(datosGA1[265:372, 2], start = 2015, frequency = 12)
IVF <- ts(datosIVF[1:108, 2], start = 2015, frequency = 12)

# Visualizacion de GA1 y IVF entre ene2015-dic2023
png("./images/EDA/Serie_GA1.png")
par(mfrow = c(1,1))
autoplot(GA1, ts.colour = "blue") + ggtitle("GA1")
dev.off()

png("./images/EDA/Serie_IVF.png")
par(mfrow = c(1,1))
autoplot(IVF, ts.colour = "black") + ggtitle("IVF")
dev.off()

# Descomposicion aditiva de las series de tiempo
png("./images/EDA/Decomp_IVF.png")
plot(decompose(IVF))
dev.off()

png("./images/EDA/Decomp_GA1.png")
plot(decompose(GA1))
dev.off()

# Graficas ACF Y PCF de series originales+
png("./images/EDA/ACF_PACF_IVF.png")
par(mfrow = c(2, 1))
acf(IVF)
pacf(IVF)
dev.off()

png("./images/EDA/ACF_PACF_GA1.png")
par(mfrow = c(2, 1))
acf(GA1)
pacf(GA1)
dev.off()
