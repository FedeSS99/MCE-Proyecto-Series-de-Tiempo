calcular_MAE_MSE <- function(Ypred, Yreal){
    DifY <- Ypred - Yreal
    MAE <- mean(abs(DifY))
    MSE <- mean(DifY^2)

    Errores <- data.frame(MAE = MAE, MSE = MSE)
    return(Errores)
}