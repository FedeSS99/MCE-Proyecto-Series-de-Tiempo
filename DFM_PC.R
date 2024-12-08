source("./cod/Config.R")
source("./cod/ETL.R")
source("./Functions_DFM.R")

# Tests aumentados de Dickey-Fuller de todas las series
adf_testGA1 <- adf.test(GA1, alternative = "stationary")
adf_testGA1
# p-valor = menor a 0.01 -> Evidencia de estacionariedad

adf_testIVF <- adf.test(IVF, alternative = "stationary")
adf_testIVF
# p-valor = 0.0153 -> Evidencia de estacionariedad

adf_testIDF1 <- adf.test(IDF1, alternative = "stationary")
adf_testIDF1
# p-valor = menor a 0.0361 -> Evidencia de estacionariedad

adf_testIDF2 <- adf.test(IDF2, alternative = "stationary")
adf_testIDF2
# p-valor = 0.1819 -> Evidencia de NO estacionariedad
adf_testDiffIDF2 <- adf.test(diff(IDF2), alternative = "stationary")
adf_testDiffIDF2
# p-valor = menor a 0.01 -> Evidencia de estacionariedad

adf_testPREC <- adf.test(PREC, alternative = "stationary")
adf_testPREC
# p-valor = menor a 0.01 -> Evidencia de estacionariedad

adf_testTEMP <- adf.test(TEMP, alternative = "stationary")
adf_testTEMP
# p-valor = 0.0153 -> Evidencia de estacionariedad

adf_testDIES <- adf.test(DIES, alternative = "stationary")
adf_testDIES
# p-valor = 0.6608 -> Evidencia de NO estacionariedad
adf_testDiffDIES <- adf.test(diff(DIES), alternative = "stationary")
adf_testDiffDIES
# p-valor = 0.05024 -> Evidencia de NO estacionariedad

# -- Creamos arreglo de regresores --
N <- length(GA1)
h <- 24
X <- cbind(IVF, IDF1, IDF2, PREC, TEMP, DIES)
Y24 <- GA1[(N - h + 1):N]

# -- Evaluamos metricas de información para R_max = P - 1 factores --
rmax <- ncol(X) - 1
par(mfrow = c(1, 1))
factors_1 <- bai_ng(X, demean = 2, k_max = rmax)
factors_1
png("./images/Resultados/ICFactores_PC.png")
ts.plot(factors_1$IC, col = 1:3, lwd = 2, xlab = "Number of factors", ylab = "Information Criteria")
legend(x = "topleft", legend = c("IC1", "IC2", "IC3"), col = 1:3, lwd = 2)
dev.off()

factors_2 <- onatski2010(X, demean = 2)
factors_2
factors_3 <- ratio.test(X, kmax = rmax, demean = 2)
factors_3

# -- Hallamos estimaciones de los factores --
num_factors <- 1
PC_Estim <- pc(X, r = num_factors, demean = 2)
PC_Estim$Phat
png("./images/Resultados/FactoresPC_DFM.png")
ts.plot(scale(cbind(GA1, PC_Estim$Fhat)), col = 1:( 1+ num_factors), lwd = 2, xlab = "Tiempo", ylab = "F1", main = "GA1 y Factor de variables regresoras")
legend(x = "topleft", legend = c("GA1", "F1"), col = 1:(1+ num_factors), lwd = 2)
dev.off()

# -- Evaluamos los supuestos --
LM_Y_factors <- lm(as.matrix(GA1) ~ PC_Estim$Fhat)
Residuales_Factores <- resid(LM_Y_factors)
png("./images/Resultados/ResidualesFactorGA1.png")
plot(Residuales_Factores, type = "b", lwd = 2, pch = 19, xlab = "Indice de Tiempo", ylab = "Magnitud de residual", main = "Residuales de F1 respecto a GA1")
dev.off()

adf.test(Residuales_Factores, k = 0)
# Contamos con un p-valor < 0.05 por lo que los residuales son estacionarios

nowcasts <- numeric(h)
for(hi in 1:h){
    PC_Estim_hi <- pc(X[1:(N + hi - h), ], r = num_factors, demean = 2)
    DFM_PC_model <- auto.arima(GA1[1:(N + hi - h - 1)], xreg = PC_Estim_hi$Fhat[1:(N + hi - h - 1), , drop = FALSE])
    DFM_PC_model

    GA1_presente_DFM_PC <- forecast(DFM_PC_model, xreg = PC_Estim_hi$Fhat[N + hi - h, , drop = FALSE], h = 1)
    nowcasts[hi] <- GA1_presente_DFM_PC$mean
    print(hi)
}

# Calculamos accuracy de los valores predichos por DFM sobre los valores reales
accuracy_PC <- accuracy(nowcasts, Y24)
accuracy_PC

png("./images/Resultados/PredReal_DFM.png")
plot(as.vector(nowcasts), Y24, type = "p", pch = 19,
main = "Comparación GA1 con modelo DFM por PC",
xlab = "GA1 predicho",
ylab = "GA1 real")
abline(0, 1, col = "red")
dev.off()

write.csv(data.frame(nowDFM = as.vector(nowcasts), Obs = Y24), file = "./data/nowDFM.csv", row.names = FALSE)

par(mar=c(5, 5, 5, 5), xpd=TRUE)
plot(GA1)
lines(ts(GA1[1:(N + 1 - h - 1)], start = c(2015), frequency = 12), col = 'blue')
lines(ts(nowcasts, frequency = 12, start = c(2022, 1)), col = 'red')
legend(x = 'topleft', inset=c(0, -0.3), legend = c("Valores reales", "Intervalo de entrenamiento compartido", "Nowcasting"),col = c('black', 'blue', 'red'), lwd = 2)
