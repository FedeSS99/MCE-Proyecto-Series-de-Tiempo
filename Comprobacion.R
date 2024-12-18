source("./FuncionesErrores.R")

ResSARIMA <- read.csv("./data/nowSARIMA.csv")
ResDFM <- read.csv("./data/nowDFM.csv")

# -- Calculo del MAE y MSE para los modelos SARIMA y DMF --
Errores_SARIMA <- calcular_MAE_RMSE(ResSARIMA$nowSARIMA, ResSARIMA$Obs)
rownames(Errores_SARIMA) <- c("SARIMA")
Errores_DFM <- calcular_MAE_RMSE(ResDFM$nowDFM, ResDFM$Obs)
rownames(Errores_DFM) <- c("DFM")
Errores <- rbind(Errores_SARIMA, Errores_DFM)
Errores

# -- Obtención de residuales por modelo para comparar visualmente y mediante el test de Diebold-Mariano --
errorSARIMA <- ResSARIMA$nowSARIMA - ResSARIMA$Obs
errorDFM <- ResDFM$nowDFM - ResDFM$Obs

png("./images/Resultados/ResidualesModelos.png")
ts.plot(cbind(errorSARIMA, errorDFM), lwd = 2, col = 1:2,
main = "Residuales de modelos SARIMA y DFM",
xlab = "Tiempo",
ylab = "Error")
legend(x = "bottomright", legend = c("SARIMA", "DFM"), col = 1:2, lwd = 2)
dev.off()

forecast::dm.test(errorDFM, errorSARIMA, alternative = "two.sided")
# p-valor = 0.0221 < 0.05 por lo que hay evidencia de que la eficiencia predictiva no es igual entre
# el modelo SARIMA y el modelo DFM
