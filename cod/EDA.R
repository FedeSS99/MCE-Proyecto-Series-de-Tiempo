#######
# EDA #
#######

# Visualizacion ene2015-dic2023
png("./images/EDA/Serie_GA1.png")
par(mfrow = c(1,1))
autoplot(GA1, ts.colour = "blue") + ggtitle("GA1")
dev.off()

png("./images/EDA/Serie_IVF.png")
par(mfrow = c(1,1))
autoplot(IVF, ts.colour = "black") + ggtitle("IVF")
dev.off()

png("./images/EDA/Serie_IDF1.png")
par(mfrow = c(1,1))
autoplot(IDF1, ts.colour = "orange") + ggtitle("CIDF en el ambiente y la gestion ambiental")
dev.off()

png("./images/EDA/Serie_IDF2.png")
par(mfrow = c(1,1))
autoplot(IDF2, ts.colour = "orange") + ggtitle("CIDF en vias de comunicación y correspondencia")
dev.off()

png("./images/EDA/Serie_PREC.png")
par(mfrow = c(1,1))
autoplot(PREC, ts.colour = "darkgreen") + ggtitle("Precipitación a nivel nacional")
dev.off()

png("./images/EDA/Serie_TEMP.png")
par(mfrow = c(1,1))
autoplot(TEMP, ts.colour = "darkgreen") + ggtitle("Temperatura a nivel nacional")
dev.off()

png("./images/EDA/Serie_DIES.png")
par(mfrow = c(1,1))
autoplot(DIES, ts.colour = "chocolate") + ggtitle("Precio del diésel")
dev.off()

# Descomposicion aditiva de las series de tiempo
png("./images/EDA/Decomp_IVF.png")
plot(decompose(IVF))
dev.off()

png("./images/EDA/Decomp_GA1.png")
plot(decompose(GA1))
dev.off()

png("./images/EDA/Decomp_IDF1.png")
plot(decompose(IDF1))
dev.off()

png("./images/EDA/Decomp_IDF2.png")
plot(decompose(IDF2))
dev.off()

png("./images/EDA/Decomp_PREC.png")
plot(decompose(PREC))
dev.off()

png("./images/EDA/Decomp_TEMP.png")
plot(decompose(TEMP))
dev.off()

png("./images/EDA/Decomp_DIES.png")
plot(decompose(DIES))
dev.off()

# Graficas ACF Y PCF de series originales
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