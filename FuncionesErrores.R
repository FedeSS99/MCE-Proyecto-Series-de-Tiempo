calcular_MAE_RMSE <- function(Ypred, Yreal){
    DifY <- Ypred - Yreal
    MAE <- mean(abs(DifY))
    RMSE <- sqrt(mean(DifY^2))

    Errores <- data.frame(MAE = MAE, RMSE = RMSE)
    return(Errores)
}