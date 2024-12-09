ResSARIMA <- read.csv("./data/nowSARIMA.csv")
ResDFM <- read.csv("./data/nowDFM.csv")

errorSARIMA <- ResSARIMA$nowSARIMA - ResSARIMA$Obs
errorDFM <- ResDFM$nowDFM - ResDFM$Obs

png("./images/Resultados/ResidualesModelos.png")
ts.plot(cbind(errorSARIMA, errorDFM), lwd = 2, col = 1:2,
main = "Residuales de modelos SARIMA y DFM",
xlab = "Tiempo",
ylab = "Error")
legend(x = "bottomright", legend = c("SARIMA", "DFM"), col = 1:2, lwd = 2)
dev.off()

forecast::dm.test(errorDFM, errorSARIMA, alternative = "greater")