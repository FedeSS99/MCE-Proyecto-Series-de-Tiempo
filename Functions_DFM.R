# Bai and Ng (2002) information criteria
# Y is the time series
# demean is the type of transformation
# k_max is the maximum number of factors to be tested
bai_ng <- function(Y, demean = 0, k_max = 8){
  # N and T
  n <- ncol(Y)
  tt <- nrow(Y)
  
  # minimum between N and T
  CNT2 <- min(c(n,tt))
  
  # information criteria
  IC <- matrix(0, k_max + 1, 3)
  colnames(IC) <- c("IC1", "IC2", "IC3")
  
  # estimating
  for(k in 0 : k_max){
    if(k > 0){
      pc_est <- pc(Y, k, demean)
      Phat <- pc_est$Phat
      Fhat <- pc_est$Fhat
      
      ehat <- Y - Fhat%*%t(Phat)
    }else{
      ehat <- Y
    }
    
    # variance of the errors
    Vfk <- mean(colMeans(ehat^2))
    
     # criteria
    IC[k+1, "IC1"] <- log(Vfk) + (k*((n+tt)/(n*tt)))*log((n*tt)/(n+tt))
    IC[k+1, "IC2"] <- log(Vfk) + (k*((n+tt)/(n*tt)))*log(CNT2)
    IC[k+1, "IC3"] <- log(Vfk) + (k*((log(CNT2)/CNT2)))
  }
      
  # r hat
  rhat <- apply(IC, 2, function(x) which(x == min(x)) - 1)
  
  # result
  result <- list(IC, rhat)
  names(result) <- c("IC", "rhat")
  
  return(result)
}

# number of factors procedure of Onatski (2010)
# X is the database
# demean is the type of transformation
onatski2010 <- function(X, demean = 2){  

  # raw data
  if(demean == 0)
    X <- X

  # demean data
  if(demean == 1)
    X <- scale(X, scale = FALSE)

  # standardized data
  if(demean == 2)
    X <- scale(X, scale = TRUE)
                    
  # transpose
  X <- t(X)
  
  # T and N
  T <- ncol(X)
  N <- nrow(X)
  
  # r max (Onatski WP)
  rmax <- round(1.55*min(T^(2/5), N^(2/5)))
    
  # j initial
  j <- rmax + 1

  # eigen value decomposition
  XXt <- X%*%t(X)
  SVD <- eigen(XXt/T)

  # lambda
  lambda <- SVD$values
  
  # inital values
  rhat <- Inf
  conver <- Inf
  i <- 1
  Conver <- Inf

  while(length(rhat) < 4 | conver > 0 | conver < 0){
    # iterative
    
    # betahat
    betahat <- 2*abs(coef(lm(lambda[j:(j+4)] ~ c(((j - 1):(j + 3))^(2/3))))[2])
    
    # criterion (page 1008)
    cond <- which(lambda[-length(lambda)] - lambda[-1] >= betahat)
            
    if(length(cond) != 0){
      rHat <- max(cond)
      rhat[i] <- rHat  

      # repeat 2 and 3
      j <- rhat[i] + 1
      # convergence                               
      if(length(rhat) > 4)
        conver <- sum(diff(rhat[(length(rhat)-4):length(rhat)]))
    }else{
      rhat <- rep(0, i)
      conver <- 0
    }    
    
    i <- i + 1
  }    

  # rhat
  r.hat <- c(rhat[length(rhat)], betahat)
  names(r.hat) <- c("ed", "betahat")
  
  # return
  return(r.hat)
}

# eigenvalue ratio test
# X is the database
# kmax is the maximum number to test
# demean is the type of transformation
ratio.test <- function(X, kmax = 8, demean = 2){

  # raw data
  if(demean == 0)
    X <- X

  # demean data
  if(demean == 1)
    X <- scale(X, scale = FALSE)

  # standardized data
  if(demean == 2)
    X <- scale(X, scale = TRUE)

  # T and N
  T <- nrow(X)
  N <- ncol(X)                      
  
  # extracting the eigenvalues
  if(T > N){
    mu <- eigen(t(X)%*%X/(N*T))$values
  }else{
    mu <- eigen(X%*%t(X)/(N*T))$values
  }
  
  # Vfk
  Vfk <- c()
  for(i in 1 : c(kmax + 1))        
    Vfk[i] <- sum(mu[i:min(N,T)])
    
  # mock eigenvalue
  mu0 <- sum(mu[1:(min(N,T))])/(log(min(N,T))*N)
  mu <- c(mu0, mu)
         
  # test ER(k)
  ER <- mu[1:c(kmax+2)][-c(kmax+2)]/mu[1:c(kmax+2)][-1]
  if(is.null(Vfk)){
    # test GR(k)
    GR <- NULL
  }else{
    # test GR(k)
    GR <- log(1 + mu[1:c(kmax+2)][-c(kmax+2)]/c(Vfk))/
          log(1 + mu[1:c(kmax+2)][-1]/c(Vfk))
  }
  # estimators
  ker <- which(ER == max(ER)) - 1
  kgr <- which(GR == max(GR)) - 1
  
  # result
  result <- c(ker, kgr, mu)
  names(result) <- c("ker", "kgr", paste("mu_",0:min(N,T), sep = ""))
  # return
  return(result)
}

# ADF test
adf <- function(x, type = c("none", "const", "trend"),
                alternative = c("stationary", "explosive"),
                k = trunc((length(x) - 1)^(1/3))){

  if (NCOL(x) > 1)
      stop("x is not a vector or univariate time series")
  if (any(is.na(x)))
      stop("NAs in x")
  if (k < 0)
      stop("k negative")

  type <- match.arg(type)
  alternative <- match.arg(alternative)
  DNAME <- deparse(substitute(x))

  k <- k + 1
  x <- as.vector(x, mode = "double")
  y <- diff(x)
  n <- length(y)
  z <- embed(y, k)
  yt <- z[, 1]
  xt1 <- x[k:n]
  tt <- k:n

  if(type == "none"){
    if(k > 1){
      yt1 <- z[, 2:k]
      res <- lm(yt ~ xt1 + yt1 - 1)
    }else{
      res <- lm(yt ~ xt1 - 1)
    }
    res.sum <- summary(res)
    STAT <- res.sum$coefficients[1, 1]/res.sum$coefficients[1,
              2]
    Table <- adfTable(trend = "nc", statistic = "t")
  }

  if(type == "const"){
    if(k > 1){
      yt1 <- z[, 2:k]
      res <- lm(yt ~ xt1 + 1 + yt1)
    }else{
      res <- lm(yt ~ xt1 + 1)
    }
    res.sum <- summary(res)
    STAT <- res.sum$coefficients[2, 1]/res.sum$coefficients[2,
              2]
    Table <- adfTable(trend = "c", statistic = "t")
  }

  if(type == "trend"){
    if (k > 1) {
      yt1 <- z[, 2:k]
      res <- lm(yt ~ xt1 + 1 + tt + yt1)
    }else{
      res <- lm(yt ~ xt1 + 1 + tt)
    }
    res.sum <- summary(res)
    STAT <- res.sum$coefficients[2, 1]/res.sum$coefficients[2,
              2]
    Table <- adfTable(trend = "ct", statistic = "t")
  }
  bic <- BIC(res)
  table <- Table$z
  tablen <- dim(table)[2]
  tableT <- Table$x
  tablep <- Table$y
  tableipl <- numeric(tablen)
  for (i in (1:tablen)) tableipl[i] <- approx(tableT, table[,
      i], n, rule = 2)$y
  interpol <- approx(tableipl, tablep, STAT, rule = 2)$y
  if (is.na(approx(tableipl, tablep, STAT, rule = 1)$y))
      if (interpol == min(tablep))
          warning("p-value smaller than printed p-value")
      else warning("p-value greater than printed p-value")
  if (alternative == "stationary")
      PVAL <- interpol
  else if (alternative == "explosive")
      PVAL <- 1 - interpol
  else stop("irregular alternative")
  PARAMETER <- k - 1
  METHOD <- "Augmented Dickey-Fuller Test"
  names(STAT) <- "Dickey-Fuller"
  names(PARAMETER) <- "Lag order"
  structure(list(bic = bic, statistic = STAT, parameter = PARAMETER,
      alternative = alternative, p.value = PVAL, method = METHOD,
      data.name = DNAME), class = "htest")
}

# ADF tables
adfTable <- function(trend = c("nc", "c", "ct"), statistic = c("t", "n"),
           includeInf = TRUE){
  # A function implemented by Diethelm Wuertz

  # Description:
  #   Tables critical values for augmented Dickey-Fuller test.

  # Note:
  #   x=-3:0; y=0:3; z=outer(x,y,"*"); rownames(z)=x; colnames(z)=y; z

  # Examples:
  #   adfTable()

  # FUNCTION:

  # Match Arguments:
  type = trend = match.arg(trend)
  statistic = match.arg(statistic)

  # Tables:
  if (statistic == "t") {
    # Hamilton Table B.6 - OLS t-Statistic
    if (type == "nc") {
      table = cbind(
        c(-2.66, -2.26, -1.95, -1.60, +0.92, +1.33, +1.70, +2.16),
        c(-2.62, -2.25, -1.95, -1.61, +0.91, +1.31, +1.66, +2.08),
        c(-2.60, -2.24, -1.95, -1.61, +0.90, +1.29, +1.64, +2.03),
        c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.29, +1.63, +2.01),
        c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00),
        c(-2.58, -2.23, -1.95, -1.62, +0.89, +1.28, +1.62, +2.00))
    } else if (type == "c") {
      table = cbind(
        c(-3.75, -3.33, -3.00, -2.63, -0.37, +0.00, +0.34, +0.72),
        c(-3.58, -3.22, -2.93, -2.60, -0.40, -0.03, +0.29, +0.66),
        c(-3.51, -3.17, -2.89, -2.58, -0.42, -0.05, +0.26, +0.63),
        c(-3.46, -3.14, -2.88, -2.57, -0.42, -0.06, +0.24, +0.62),
        c(-3.44, -3.13, -2.87, -2.57, -0.43, -0.07, +0.24, +0.61),
        c(-3.43, -3.12, -2.86, -2.57, -0.44, -0.07, +0.23, +0.60))
    } else if (type == "ct") {
      table = cbind(
        c(-4.38, -3.95, -3.60, -3.24, -1.14, -0.80, -0.50, -0.15),
        c(-4.15, -3.80, -3.50, -3.18, -1.19, -0.87, -0.58, -0.24),
        c(-4.04, -3.73, -3.45, -3.15, -1.22, -0.90, -0.62, -0.28),
        c(-3.99, -3.69, -3.43, -3.13, -1.23, -0.92, -0.64, -0.31),
        c(-3.98, -3.68, -3.42, -3.13, -1.24, -0.93, -0.65, -0.32),
        c(-3.96, -3.66, -3.41, -3.12, -1.25, -0.94, -0.66, -0.33))
    } else {
      stop("Invalid type specified")
    }
  } else if (statistic == "z" || statistic == "n") {
    # Hamilton Table B.5 - Based on OLS Autoregressive Coefficient
    if (type == "nc") {
      table = cbind(
        c(-11.9,  -9.3,  -7.3,  -5.3, +1.01, +1.40, +1.79, +2.28),
        c(-12.9,  -9.9,  -7.7,  -5.5, +0.97, +1.35, +1.70, +2.16),
        c(-13.3, -10.2,  -7.9,  -5.6, +0.95, +1.31, +1.65, +2.09),
        c(-13.6, -10.3,  -8.0,  -5.7, +0.93, +1.28, +1.62, +2.04),
        c(-13.7, -10.4,  -8.0,  -5.7, +0.93, +1.28, +1.61, +2.04),
        c(-13.8, -10.5,  -8.1,  -5.7, +0.93, +1.28, +1.60, +2.03))
    } else if (type == "c") {
      table = cbind(
        c(-17.2, -14.6, -12.5, -10.2, -0.76, +0.01, +0.65, +1.40),
        c(-18.9, -15.7, -13.3, -10.7, -0.81, -0.07, +0.53, +1.22),
        c(-19.8, -16.3, -13.7, -11.0, -0.83, -0.10, +0.47, +1.14),
        c(-20.3, -16.6, -14.0, -11.2, -0.84, -0.12, +0.43, +1.09),
        c(-20.5, -16.8, -14.0, -11.2, -0.84, -0.13, +0.42, +1.06),
        c(-20.7, -16.9, -14.1, -11.3, -0.85, -0.13, +0.41, +1.04))
    } else if (type == "ct") {
      table = cbind(
        c(-22.5, -19.9, -17.9, -15.6, -3.66, -2.51, -1.53, -0.43),
        c(-25.7, -22.4, -19.8, -16.8, -3.71, -2.60, -1.66, -0.65),
        c(-27.4, -23.6, -20.7, -17.5, -3.74, -2.62, -1.73, -0.75),
        c(-28.4, -24.4, -21.3, -18.0, -3.75, -2.64, -1.78, -0.82),
        c(-28.9, -24.8, -21.5, -18.1, -3.76, -2.65, -1.78, -0.84),
        c(-29.5, -25.1, -21.8, -18.3, -3.77, -2.66, -1.79, -0.87))
    } else {
      stop("Invalid type specified")
    }
  } else {
    stop("Invalid statistic specified")
  }

  # Transpose:
  Table = t(table)
  colnames(Table) = c("0.010", "0.025", "0.050", "0.100", "0.900",
                      "0.950", "0.975", "0.990")
  rownames(Table) = c(" 25", " 50", "100", "250", "500", "Inf")
  ans = list(
    x = as.numeric(rownames(Table)),
    y = as.numeric(colnames(Table)),
    z = Table)
  class(ans) = "gridData"

  # Exclude Inf:
  if (!includeInf) {
    nX = length(ans$x)
    ans$x = ans$x[-nX]
    ans$z = ans$z[-nX, ]
  }

  # Add Control:
  attr(ans, "control") <-
    c(table = "adf", trend = trend, statistic = statistic)

  # Return Value:
  ans
}

# principal components factor extraction
# Y is the matrix of time series
# r is the number of factors
# demean is the type of transformation
pc <- function(Y, r = 1, demean = 2){  
  # number of time series and observations
  N <- ncol(Y)
  Ty <- nrow(Y)
  
  # transformation
  if(demean == 0)
    Y <- Y
  
  if(demean == 1)
    Y <- scale(Y, center = TRUE)
  
  if(demean == 2)
    Y <- scale(Y)
  
  # eigenvalue decomposition
  ed <- eigen(t(Y)%*%Y)
  
  # loading matrix
  Phat <- sqrt(N)*ed$vectors[,1:r,drop=FALSE] 
  rownames(Phat) <- colnames(Y)
  
  # factor estimates
  Fhat <- Y%*%Phat/N
  
  # idiosyncratic errors
  ehat <- Y - Fhat%*%t(Phat)
  
  colnames(Phat) <- colnames(Fhat) <- paste("Factor",1:r, sep = "")
  
  # result
  result <- list(ed, Y, Phat, Fhat, ehat)
  names(result) <- c("ed", "Y", "Phat", "Fhat", "ehat")

  return(result)
}

# confidence intervals for pc factor extraction
# pc_obj is the pc factor extraction
# x is a target time series
# alpha is the signficance level
# qs parameter of autocorrelation
pc_ci <- function(pc_obj, x = NULL, alpha = 0.05, qs = 2){
  # objects
  ed <- pc_obj$ed
  Y <- pc_obj$Y
  Phat <- pc_obj$Phat
  Fhat <- pc_obj$Fhat
  ehat <- pc_obj$ehat
  
  r <- ncol(Fhat)
  N <- ncol(Y)
  Ty <- nrow(Y)
  
  # rotating according to target time series  
  if(!is.null(x)){
    rho <- cor(x, Fhat, use = "pairwise.complete.obs")
  
    Phat[,rho < 0] <- -Phat[,rho < 0]
    Fhat[,rho < 0] <- -Fhat[,rho < 0]
  }
  
  # estimating confidence intervals
  V <- diag(ed$values[1:r],r,r)
  Vinv <- solve(V/(N*Ty))
  
  # z value
  z <- qnorm(1-alpha/2, 0, 1)
  
  # @ static factor confidence interval 
  cil <- ciu <- St <- matrix(0, Ty, r)
  
  # asymptotic variance for factors
  for (t in 1:Ty) {
    Gamma <- matrix(0, r, r)
  
    for (i in 1:N) {
      Gamma.i <- matrix(Phat[i,])%*%t(matrix(Phat[i,]))*(ehat[t,i]^2)
      Gamma <- Gamma + Gamma.i
    }
    Gamma <- (1/N)*(Gamma)
    St[t,] <- sqrt(Vinv%*%Gamma%*%Vinv)
  
    # confidence Intervals
    cil[t,] <- Fhat[t,] - z*N^(-1/2)*St[t,]
    ciu[t,] <- Fhat[t,] + z*N^(-1/2)*St[t,]
  }
  

  # asymptotic variance for loadings
  Phi_list <- list()
  for(i in 1 : N){
    Phi <- matrix(0, r, r)
    for (t in 1:Ty) {
      Phi.i <- matrix(Fhat[t,])%*%t(matrix(Fhat[t,]))*(ehat[t,i]^2)
      Phi <- Phi + Phi.i
    }
    Phi <- Phi/Ty
    Phi_list[[i]] <- Phi
  }
  
  Divq <- list()
  for(v in 1 : qs){
    Div_list <- list()
    for(i in 1 : N){
      Div <- matrix(0, r, r)
      for(t in (v+1):Ty){
        Div.i <- Fhat[t,]*Fhat[t-v,]*ehat[t,i]*ehat[t-v,i]
        Div <- Div + Div.i
      }
      Div <- (1-v/(qs+1))*(Div + t(Div))
      Div_list[[i]] <- Div/Ty
    }
    Divq[[v]] <- Div_list
  }
  
  D_list <- list()
  for(i in 1 : N)
    D_list[[i]] <- Divq[[1]][[i]] + Divq[[2]][[i]]
  
  se <- matrix(0, N, r)
  for(i in 1 : N){
    se[i,] <- sqrt(diag((Phi_list[[i]] + D_list[[i]])/Ty))
  }
  
  # matrices with intervals
  F_mat <- cbind(cil, Fhat, ciu)
  colnames(F_mat) <- c("lower", "Fhat", "upper")
  
  li <- Phat - z*se
  ls <- Phat + z*se
  
  P_mat <- cbind(li, Phat, ls)
  colnames(P_mat) <- c("lower", "Phat", "upper")
 
  # return
  result <- list(P_mat, F_mat)
  names(result) <- c("P_mat", "F_mat")
  
  return(result)
}


# @ Kalman smoothing factor extraction
# Y_t = PF_t + epsilon_t
# F_t = PhiF[t-1] + eta_t
# Y is the time series
# r is the number of common factors
# demean is the transformation type
# type_adf is the specification in ADF test
# lag_max is the maximum lags to test
KS <- function(Y, r, demean = 2, type_adf = "const", lag_max = 8,
               alpha = 0.05){
    
  # type of demean data 
  if(demean == 0)
    Y <- Y
  
  if(demean == 1)
    Y <- scale(Y, center = TRUE, scale = FALSE)
  
  if(demean == 2)
    Y <- scale(Y)
  
  # number of rows
  tt <- nrow(Y)
  
  # static PC factor extraction
  pc_est <- pc(Y, r, demean)
  
  # loadings and factor estimates
  Phat <- pc_est$Phat
  Fhat <- pc_est$Fhat
  
  # selecting the lags in ADF test
  BIC_factors <- matrix(0, lag_max + 1, r)
  for(i in 1 : r){
    for(j in 0 : lag_max){
      BIC_factors[j+1,i] <- adf(Fhat[,i], type_adf, k = j)$bic
    }
  }
  
  # lags
  lags <- apply(BIC_factors, 2, function(x) which(x == min(x)) - 1)
  
  # ADF tests
  adf_tests <- sapply(1:r, function(x) adf(Fhat[,x], type_adf, 
                                           k = lags[x])$p.value)
  
  # which factors are non-stationary?
  ns_factors <- adf_tests > alpha
  
  # Kalman filter parameters
  Phi <- Sigma_eta <- Sigma_eta_0 <- matrix(0, r, r)
  eta <- matrix(0, tt-1, r)
  
  # for non-stationary common factors 
  if(sum(ns_factors) > 0){
    for(i in which(ns_factors)){
      Phi[i,i] <- 1
      eta[,i] <- diff(Fhat[,i])
      Sigma_eta[i,i] <- var(eta[,i])
      Sigma_eta_0[i,i] <- 10^7 
    }
  }
  
  # for stationary common factors
  if(sum(!ns_factors) > 0){
    for(i in which(!ns_factors)){
      ar_model <- ar.ols(Fhat[,i], FALSE, 1, demean = FALSE) 
      Phi[i,i] <- ar_model$ar[,,1]
      eta[,i] <- matrix(ar_model$resid[-1])
      Sigma_eta[i,i] <- var(eta[,i])
      Sigma_eta_0[i,i] <- Sigma_eta[i,i]/(1-Phi[i,i]^2)
    }
  }
  
  # if variance is negative by Phi tending to 1
  Sigma_eta_0[Sigma_eta_0 < 0] <- 10^7
  
  # Sigma_epsilon
  chi <- Fhat%*%t(Phat)
  ehat <- Y - chi
  Sigma_epsilon <- diag(diag(cov(ehat)))
  
  # initial conditions in factors
  F0 <- colMeans(Fhat)
  
  # state-space representation
  kf <- dlm(FF = Phat, V = Sigma_epsilon, GG = Phi, W = Sigma_eta, 
            m0 = F0, C0 = Sigma_eta_0)
  
  # smoothing factors
  z <- qnorm(1-alpha/2, 0, 1)
  smooth <- dlmSmooth(Y, kf)
  
  if(r > 1){
    Fhat_s <- as.matrix(smooth$s[-1,1:r,drop=FALSE])
    Lower_s <- Fhat_s - z*as.matrix(smooth$D.S[-1,1:r,drop=FALSE]) 
    Upper_s <- Fhat_s + z*as.matrix(smooth$D.S[-1,1:r,drop=FALSE])
  }else{
    Fhat_s <- as.matrix(smooth$s[-1])
    Lower_s <- Fhat_s - z*smooth$D.S[-1]
    Upper_s <- Fhat_s + z*smooth$D.S[-1]
  }
  colnames(Fhat_s) <- colnames(Lower_s) <- colnames(Upper_s) <- colnames(Fhat)
  
  # return
  result <- list(Phat, Fhat, Fhat_s, Lower_s, Upper_s)
  names(result) <- c("Phat", "Fhat", "Fhat_s", "Lower_s", "Upper_s")
  
  return(result)
}

# @ Barigozzi et al. (2016) method
# Y is the time series
# r is the number of common factors
# demean is the transformation type
barigozzi <- function(Y, r, demean = 2){
  # number of time series
  N <- ncol(Y)
    
  # transformation
  if(demean == 0)
    Y <- Y
  
  if(demean == 1)
    Y <- scale(Y, center = TRUE)
  
  if(demean == 2)
    Y <- scale(Y)

  # principal components factor extraction with first differenced data
  pc_est <- pc(diff(Y), r = r, demean = 0)
  
  # loadings estimates
  Phat <- pc_est$Phat
  
  # factor estimates
  Fhat <- Y%*%Phat/N
  
  rownames(Phat) <- colnames(Y)
  colnames(Phat) <- colnames(Fhat)
  rownames(Fhat) <- rownames(Y)
  
  # result
  result <- list(Phat, Fhat)
  names(result) <- c("Phat", "Fhat")
  
  # return
  return(result)
}

# multivariate r2
# F is the simulated common factors
# Fhat is the estimated common factors
Tr <- function(F, Fhat){
        sum(diag(t(F)%*%Fhat%*%solve(t(Fhat)%*%Fhat)%*%t(Fhat)%*%F))/
        sum(diag(t(F)%*%F))
}

# percentage variations
# x is the time series
# q is the parameter of % variation: 1 is monthly, 12 is annual
cpm <- function(x, q = 1){

  if(q != 0){
    index <- matrix(0, nrow(x) - q, 2)
    index[,1] <- 1:nrow(index)
    index[,2] <- (q+1):(nrow(index)+q)
  
    xm <- matrix(0, nrow(index), 1)
    for(i in 1 : nrow(xm))
      xm[i,] <- x[index[i,2],]/x[index[i,1],]*100-100
    
    rownames(xm) <- rownames(x)[-(1:q)]
    colnames(xm) <- colnames(x)
  }else{
    xm <- x
  }
  return(xm)
}

# lag function
# x is the series
# lag is the lag parameter
my_lag <- function(x, lag = 1){
  if(lag != 0){
    tt <- nrow(x)
    xl <- c(rep(NA, lag), x[1:(tt-lag),])
    xl <- as.matrix(xl)
    colnames(xl) <- colnames(x)
    rownames(xl) <- rownames(x)
  }else{
    xl <- x
  }  
  return(xl)
}


# optimal transformation
# x is the target time series
# y is the y_i time series of the DFM
# start is the starting date
transformation <- function(x, y, start = "2004/01"){

  # we consider three possible transformations: monthly, annual or levels  
  y_mv <- cpm(y, 1)
  y_av <- cpm(y, 12)
  y_levels <- cpm(y, 0)
  y_lag <- my_lag(y, 1)
  
  # dates indicator
  ind_dates <- rownames(x)[which(rownames(x) == start):nrow(x)]
  
  # estimating the linear correlations
  rho_mv <- cor(x[ind_dates,,drop=FALSE], y_mv[ind_dates,,drop=FALSE],
              use = "pairwise.complete.obs")
  
  rho_av <- cor(x[ind_dates,,drop=FALSE], y_av[ind_dates,,drop=FALSE],
              use = "pairwise.complete.obs")
  
  rho_levels <- cor(x[ind_dates,,drop=FALSE], y_levels[ind_dates,,drop=FALSE],
                  use = "pairwise.complete.obs")
  
  rho_lag <- cor(x[ind_dates,,drop=FALSE], y_lag[ind_dates,,drop=FALSE],
                  use = "pairwise.complete.obs")
        
  # obtaining the maximum value
  if(abs(rho_mv) > abs(rho_av) & abs(rho_mv) > abs(rho_levels) & 
     abs(rho_mv) > abs(rho_lag))
    y_trans <- y_mv
  
  if(abs(rho_av) > abs(rho_mv) & abs(rho_av) > abs(rho_levels) & 
     abs(rho_av) > abs(rho_lag))
    y_trans <- y_av
  
  if(abs(rho_levels) > abs(rho_mv) & abs(rho_levels) > abs(rho_av) & 
     abs(rho_levels) > abs(rho_lag))
    y_trans <- y_levels
  
  if(abs(rho_lag) > abs(rho_mv) & abs(rho_lag) > abs(rho_av) & 
     abs(rho_lag) > abs(rho_levels))
    y_trans <- y_lag
  
  # return
  rho <- c(rho_levels, rho_mv, rho_av, rho_lag)
  names(rho) <- c("levels", "mv", "av", "lag")
  
  result <- list(rho, y_trans)
  names(result) <- c("rho", "trans")
  
  return(result)
}

# selecting Google Trends topics
# x is the target time series
# GT is the Google Trends database
# X is the DFM time series
# Ng is the number of topics to be tested
# size is the percentage of sample in the rolling window
# restric_X is whether you consider the possibility to non include Google Topics
pls_selection <- function(x, GT = google, X = dfm_t, Ng = 4, size = 0.7,
                          restrict_X = TRUE){

  # as matrix
  GT <- as.matrix(GT)
  
  # target time series
  xt <- x[!is.na(x),,drop=FALSE]
  
  GT <- GT[1:nrow(xt),]
  X <- X[1:nrow(xt),]
  
  # explanation variance
  var_x <- c(plsr(xt ~ X)$Xvar/plsr(xt ~ X)$Xtotvar)[1]
  
  # sample size
  sample <- round(nrow(xt)*size)
  
  # sequence of rolling sample
  S <- nrow(xt) - sample
  
  # loading matrix
  loadings_pls <- matrix(0, S, ncol(GT))
  colnames(loadings_pls) <- colnames(GT)
  
  # rolling sample PLS
  for(t in 1 : S){ 
    regre_pls_t <- plsr(xt[t:(sample+t),] ~ GT[t:(sample+t),])
    loadings_pls[t,] <- loadings(regre_pls_t)[,1]
  }
  
  # 95% confidence interval
  coef_pls <- cbind(apply(loadings_pls, 2, function(x) quantile(x, 0.025)),
                    apply(loadings_pls, 2, function(x) quantile(x, 0.5)),
                    apply(loadings_pls, 2, function(x) quantile(x, 0.975)))
  
  colnames(coef_pls) <- c("Lower", "Median", "Upper")
  
  
  coef_pls <- coef_pls[order(abs(coef_pls[,"Median"]), decreasing = TRUE),]
  
  
  # significant Google Topics by PLS
  significance <- rownames(coef_pls)[apply(coef_pls, 1, 
                          function(x) all(x < 0) | all(x > 0))]
  
  # restricting by number of topics
  if(!is.null(Ng)){ 
    if(Ng < length(significance))
      significance <- significance[1:Ng]
  }  
  
  # determining the combination of topics that maximize the sample effective 
  # dependence
  n <- length(significance)
  
  # all combinations
  var_list <- list()
  
  # searching topics
  for(i in 1 : n){ 
    Combn <- combn(n, i)
    var_vector <- rep(0, ncol(Combn))
    for(j in 1 : length(var_vector)){ 
      Y <- cbind(X, GT[,significance[Combn[,j]]])
      var_vector[j] <- c(plsr(xt ~ Y)$Xvar/plsr(xt ~ Y)$Xtotvar)[1]
    }
    var_list[[i]] <- var_vector 
  }
  
  
  # maximum variance
  max_var <- max(unlist(var_list))
  
  # number of topics
  ni <- which(sapply(1:n, 
                     function(x) max(var_list[[x]])) == max_var)
  
  # selecting
  topics <- significance[combn(n, ni)[,var_list[[ni]] == max_var]]
  
  if(restrict_X){
    if(max_var < var_x)
      topics <- "None"
  }
  
  # return
  result <- list(coef_pls, topics)
  names(result) <- c("coef_pls", "topics")
  
  return(result)
}

# introducing rapid estimates
gdp_o <- function(x_levels, target, gdp, gdp_levels, H, q, dates){

  # length
  Tn <- nrow(x_levels) - H 

  # restriction
  month_nowcast <- substring(rownames(x_levels)[Tn+1], 6, 7)
  if(month_nowcast == "03" | month_nowcast == "06" |
     month_nowcast == "09" | month_nowcast == "12"){
    
    # annual growth timely estimates
    if(month_nowcast == "03")
      quarter <- "01"

    if(month_nowcast == "06")
      quarter <- "02"

    if(month_nowcast == "09")
      quarter <- "03"

    if(month_nowcast == "12")
      quarter <- "04"
    
    # identifying the quarters
    year_nowcast <- substring(rownames(x_levels)[Tn+1], 1, 4) 
    
    current_quarter <- paste(year_nowcast,"/", quarter, sep = "") 
    previous_quarter <- paste(as.numeric(year_nowcast)-1,"/", 
                              quarter, sep = "")  
        
    # current growth
    growth <- gdp[current_quarter,target,drop=FALSE]
    
    # recover the levels
    gdp_nowcast <- gdp_levels[previous_quarter,target,drop=FALSE]*
                    (1 + growth/100)
    rownames(gdp_nowcast) <- rownames(growth)
    
    gdp_target <- rbind(gdp_levels[,target,drop=FALSE], gdp_nowcast)
    
    # previous quarterly levels
    qp <- x_levels[(Tn-11):(Tn-13),target,drop=FALSE]
    
    # current quarterly - 1 levels
    qc <- x_levels[Tn:(Tn-1),target,drop=FALSE]
    
    # Corona et al. (2024) page 12
    x_nowcast <- (1+growth/100)*sum(qp)-sum(qc)
    rownames(x_nowcast) <- rownames(x_levels)[Tn+1]
    
    # recover the levels 
    x_levels[rownames(x_nowcast),] <- x_nowcast
    
    # apply transformation
    x_trans <- cpm(x_levels, q)[dates,,drop=FALSE]
    colnames(x_trans) <- paste(target, "_O", sep = "")
  
    x_trans <- na_locf(x_trans)
  }else{
    x_trans <- NULL
  }
  return(x_trans)
}

# procedure to estimate a correlated factor with respect to a target time series
# x is the target time series
# Y is the database of DFM
# init is the parameter to starting to estimate the linear correlations
new_factor <- function(x, Y, init = 12, expansive = TRUE){

  # estimating the recursive correlations
  Rho <- matrix(NA, nrow(x)-init, ncol(Y))
  colnames(Rho) <- colnames(Y)
  
  # expansive window
  if(expansive){
    for(j in 1 : ncol(Y)){
      rho <- c()
      for(i in 1 : c(nrow(x)-init)){ 
        rho[i] <- cor(x[(nrow(x)-init+1-i):nrow(x),,drop=FALSE],
                      Y[(nrow(x)-init+1-i):nrow(x),j], 
                      use = "pairwise.complete.obs") 
      }
      Rho[,j] <- rev(rho)
    }
  }else{ # rolling window
    for(j in 1 : ncol(Y)){
      rho <- c()
      for(i in 0 : c(nrow(x)-init)){
        rho[i] <- cor(x[(nrow(x)-init+1-i):(nrow(x)-i),], 
                      Y[(nrow(x)-init+1-i):(nrow(x)-i), j], 
                      use = "pairwise.complete.obs")
      }
        Rho[,j] <- rev(rho)
    }
  }
  # completing the correlations
  Rho_init <- matrix(NA, init, ncol(Y))
  for(i in 1 : init)
    Rho_init[i,] <- Rho[1,]
  
  Rho_f <- rbind(Rho_init, Rho)
  rownames(Rho_f) <- rownames(Y)
  
  # estimating the latent time series
  f <- scale(as.matrix(rowSums(sapply(1:ncol(Y), 
                    function(x) Rho_f[,x]*scale(Y)[,x]), na.rm = TRUE)))
  rownames(f) <- rownames(Y)
  colnames(f) <- "f"
  
  # result
  result <- list(Rho_f, f)
  names(result) <- c("Rho", "f")
  
  return(result)
}

# @ different transformations of nowcasting
# ioae is the nowcasts
# target is the name of the time series
# trans is the type of transformation
# H is the nowcasting horizon
nowcasting_trans <- function(ioae, yfcst, target, trans, H){
  if(trans == "MV"){
    ioae_mv <- ioae
    
    # obtaining the levels and annual variations
    ioae_levels <- ioae_av <- matrix(NA, H, 3)
    colnames(ioae_levels) <- colnames(ioae_av) <- colnames(ioae_mv)
    rownames(ioae_levels) <- rownames(ioae_av) <- rownames(ioae_mv)
    
    # recovering the mean
    levels <- yfcst[,target,drop=FALSE]
    
    for(h in 1 : H){  
      levels[nrow(levels)-H+h,] <- levels[nrow(levels)-H+h-1,]*
                                      (1+ioae_mv[h, "mean"]/100)
      ioae_levels[h, "mean"] <- levels[nrow(levels)-H+h,]
    }
    
    # recovering the lower and upper intervals
    for(h in 1 : H){
      ioae_levels[h, "lower"] <- levels[nrow(levels)-H+h-1,]*
                                        (1+ioae_mv[h, "lower"]/100)
      
      ioae_levels[h, "upper"] <- levels[nrow(levels)-H+h-1,]*
                                        (1+ioae_mv[h, "upper"]/100)

    }
    
    # @ annual variations
    ioae_av[, "mean"] <- cpm(levels, 12)[rownames(ioae_mv),]
    
    # lower and upper intervals
    for(h in 1 : H){ 
      ioae_av[h, "lower"] <- (ioae_levels[h,"lower"]/
                              levels[nrow(levels)-H+h-12,])*100-100
  
      ioae_av[h, "upper"] <- (ioae_levels[h,"upper"]/
                              levels[nrow(levels)-H+h-12,])*100-100
    }
    
  }
  
  
  if(trans == "AV"){
    ioae_av <- ioae
    
    # obtaining the levels and monthly variations
    ioae_levels <- ioae_mv <- matrix(NA, H, 3)
    colnames(ioae_levels) <- colnames(ioae_mv) <- colnames(ioae_av)
    rownames(ioae_levels) <- rownames(ioae_mv) <- rownames(ioae_av)
    
    # recovering the mean
    levels <- yfcst[,target,drop=FALSE]
    
    for(h in 1 : H){
      levels[nrow(levels)-H+h,] <- levels[nrow(levels)-H+h-12,]*
                                      (1+ioae_av[h, "mean"]/100)
      ioae_levels[h, "mean"] <- levels[nrow(levels)-H+h,]
    }
    
    # recovering the lower and upper intervals
    for(h in 1 : H){
      ioae_levels[h, "lower"] <- levels[nrow(levels)-H+h-12,]*
                                        (1+ioae_av[h, "lower"]/100)
      
      ioae_levels[h, "upper"] <- levels[nrow(levels)-H+h-12,]*
                                        (1+ioae_av[h, "upper"]/100)

    }
    
    # @ monthly variations
    ioae_mv[, "mean"] <- cpm(levels, 1)[rownames(ioae_av),]
    
    # lower and upper intervals
    for(h in 1 : H){ 
      ioae_mv[h, "lower"] <- (ioae_levels[h,"lower"]/
                              levels[nrow(levels)-H+h-1,])*100-100
  
      ioae_mv[h, "upper"] <- (ioae_levels[h,"upper"]/
                              levels[nrow(levels)-H+h-1,])*100-100
    }
  }
  
  if(trans == "LEVELS"){
    ioae_levels <- ioae
    
    # obtaining the annual and monthly variations
    ioae_av <- ioae_mv <- matrix(NA, H, 3)
    colnames(ioae_av) <- colnames(ioae_mv) <- colnames(ioae_levels)
    rownames(ioae_av) <- rownames(ioae_mv) <- rownames(ioae_levels)
    
    # @ monthly variations
    levels <- yfcst[,target,drop=FALSE]

    levels[is.na(levels)] <- ioae_levels[,"mean"]
    ioae_mv[, "mean"] <- cpm(levels, 1)[rownames(ioae_levels),]
    
    # lower and upper intervals
    Ty_h <- nrow(levels) - H
    
    for(h in 1 : H){
      ioae_mv[h, "lower"] <- (ioae_levels[h,"lower"]/
                              levels[nrow(levels)-H+h-1,])*100-100
  
      ioae_mv[h, "upper"] <- (ioae_levels[h,"upper"]/
                              levels[nrow(levels)-H+h-1,])*100-100
    }

    
    # @ annual variations
    levels <- yfcst[,target,drop=FALSE]

    levels[is.na(levels)] <- ioae_levels[,"mean"]
    ioae_av[, "mean"] <- cpm(levels, 12)[rownames(ioae_levels),]
    
    for(h in 1 : H){
      ioae_av[h, "lower"] <- (ioae_levels[h,"lower"]/
                              levels[nrow(levels)-H+h-12,])*100-100
  
      ioae_av[h, "upper"] <- (ioae_levels[h,"upper"]/
                              levels[nrow(levels)-H+h-12,])*100-100
    }

  }
  
  # result
  result <- list(ioae_mv, ioae_av, ioae_levels)
  names(result) <- c("MV", "AV", "LEVELS")
  
  return(result)
}


# trans is the type of transformation
# Ht is the sample size of CV
# H is the nowcast horizon
# x is the target time series
# out_sample is the indicator of out sample CV
# yfcst is the matrix of target time series
# list_mean includes the historical nowcasts
# list_lower includes the historical lower intervals
# list_upper includes the historical upper intervals
# mae_list is the object of MAE
# MAE_h is the Mean Absolute Errors
# MedAE is the median of Absolute Errors
alternative_nowcasts <- function(trans, Ht, H, x, out_sample, yfcst,
                                 list_mean, list_lower, list_upper,
                                 mae_list, MAE_h, MedAE){
  # if monthly variations
  if(trans == "MV"){
    # annual variations
    mae_av_list <- list()
    for(j in 1 : H){ 
        
      # matrix of annual variations
      mae_av_matrix <- matrix(NA, Ht, 4)
      colnames(mae_av_matrix) <- c("observed", "lower", "mean", "upper")
      rownames(mae_av_matrix) <- rownames(xt)[out_sample[,j]]
      
      # results along the time
      for(h in 1 : Ht){ 
        
        # months to predict
        months_h <- rownames(x)[out_sample[h,]]
      
        # levels
        yfcst_h <- yfcst[1:which(rownames(yfcst) == months_h[H]),,drop=FALSE]
      
        rownames(list_mean[[h]]) <- rownames(list_lower[[h]]) <- 
          rownames(list_upper[[h]]) <- months_h
      
        levels_h <- yfcst[months_h, target, drop=FALSE]
        
        # annual variations
        av_h <- cpm(yfcst_h[,target,drop=FALSE], 12)[months_h,,drop=FALSE]
    
        # nowcasts
        ioae_h <- cbind(list_lower[[h]][,j], list_mean[[h]][,j],
                          list_upper[[h]][,j])
    
        colnames(ioae_h) <- c("lower", "mean", "upper")
        rownames(ioae_h) <- months_h
        
        Ty_h <- nrow(yfcst_h) - H
        yfcst_h[(Ty_h+1):(Ty_h+H), ] <- NA
        
        # transformations
        nowcast_h <- nowcasting_trans(ioae_h,
                            yfcst_h, target, trans, H)
      
        # final result
        mae_av_matrix[h, ] <- c(av_h[j,], nowcast_h$AV[j,])
      }
      mae_av_list[[j]] <- mae_av_matrix
    }
    # statistics
    MAE_av_h <- colMeans(sapply(1:H, 
                                function(x) abs(mae_av_list[[x]][,"observed"] -
                                  mae_av_list[[x]][,"mean"])))
    MedAE_av_h <- apply(sapply(1:H, 
                                function(x) abs(mae_av_list[[x]][,"observed"] -
                                  mae_av_list[[x]][,"mean"])), 2, median)
    # levels
    mae_levels_list <- list()
    for(j in 1 : H){ 
        
      # see previous comments (first example)
      mae_levels_matrix <- matrix(NA, Ht, 4)
      colnames(mae_levels_matrix) <- c("observed", "lower", "mean", "upper")
      rownames(mae_levels_matrix) <- rownames(xt)[out_sample[,j]]
      
      for(h in 1 : Ht){
      
        # see previous comments (first example)  
        months_h <- rownames(x)[out_sample[h,]]
      
        yfcst_h <- yfcst[1:which(rownames(yfcst) == months_h[H]),,drop=FALSE]
      
        rownames(list_mean[[h]]) <- rownames(list_lower[[h]]) <- 
          rownames(list_upper[[h]]) <- months_h
      
        # see previous comments (first example)
        levels_h <- yfcst[months_h, target, drop=FALSE]
    
        ioae_h <- cbind(list_lower[[h]][,j], list_mean[[h]][,j],
                          list_upper[[h]][,j])
    
        colnames(ioae_h) <- c("lower", "mean", "upper")
        rownames(ioae_h) <- months_h
    
        Ty_h <- nrow(yfcst_h) - H
        yfcst_h[(Ty_h+1):(Ty_h+H), ] <- NA
        
        # see previous comments (first example)
        nowcast_h <- nowcasting_trans(ioae_h,
                            yfcst_h, target, trans, H)
      
        mae_levels_matrix[h, ] <- c(levels_h[j,], nowcast_h$LEVELS[j,])
      }
      mae_levels_list[[j]] <- mae_levels_matrix
    }
    
    # see previous comments (first example)
    MAE_levels_h <- colMeans(sapply(1:H, 
                                function(x) abs(mae_levels_list[[x]][,
                                "observed"] - mae_levels_list[[x]][,"mean"])))
    
    MedAE_levels_h <- apply(sapply(1:H, 
                                function(x) abs(mae_levels_list[[x]][,
                                "observed"] - mae_levels_list[[x]][,"mean"])),
                                2, median)
    
    # original tables
    mae_mv_list <- mae_list  
    MAE_mv_h <- MAE_h
    MedAE_mv_h <- MedAE
  }
  
  if(trans == "AV"){
    # monthly variations
    mae_mv_list <- list()
    for(j in 1 : H){ 
        
      mae_mv_matrix <- matrix(NA, Ht, 4)
      colnames(mae_mv_matrix) <- c("observed", "lower", "mean", "upper")
      rownames(mae_mv_matrix) <- rownames(xt)[out_sample[,j]]
      
      for(h in 1 : Ht){ 
        # see previous comments (first example)
        months_h <- rownames(x)[out_sample[h,]]
      
        yfcst_h <- yfcst[1:which(rownames(yfcst) == months_h[H]),,drop=FALSE]
      
        rownames(list_mean[[h]]) <- rownames(list_lower[[h]]) <- 
          rownames(list_upper[[h]]) <- months_h
      
        # see previous comments (first example)
        levels_h <- yfcst[months_h, target, drop=FALSE]
        mv_h <- cpm(yfcst_h[,target,drop=FALSE], 1)[months_h,,drop=FALSE]
    
      
        ioae_h <- cbind(list_lower[[h]][,j], list_mean[[h]][,j],
                          list_upper[[h]][,j])
    
        colnames(ioae_h) <- c("lower", "mean", "upper")
        rownames(ioae_h) <- months_h
    
        Ty_h <- nrow(yfcst_h) - H
        yfcst_h[(Ty_h+1):(Ty_h+H), ] <- NA

        # see previous comments (first example)
        nowcast_h <- nowcasting_trans(ioae_h,
                            yfcst_h, target, trans, H)
      
        mae_mv_matrix[h, ] <- c(mv_h[j,], nowcast_h$MV[j,])
      }
      mae_mv_list[[j]] <- mae_mv_matrix
    }
    
    # see previous comments (first example)
    MAE_mv_h <- colMeans(sapply(1:H, 
                                function(x) abs(mae_mv_list[[x]][,"observed"] -
                                  mae_mv_list[[x]][,"mean"])))
    MedAE_mv_h <- apply(sapply(1:H, 
                                function(x) abs(mae_mv_list[[x]][,"observed"] -
                                  mae_mv_list[[x]][,"mean"])), 2, median)
    
    # levels
    mae_levels_list <- list()
    for(j in 1 : H){ 
        
      # see previous comments (first example)
      mae_levels_matrix <- matrix(NA, Ht, 4)
      colnames(mae_levels_matrix) <- c("observed", "lower", "mean", "upper")
      rownames(mae_levels_matrix) <- rownames(xt)[out_sample[,j]]
      
      for(h in 1 : Ht){
        # see previous comments (first example)
        months_h <- rownames(x)[out_sample[h,]]
      
        yfcst_h <- yfcst[1:which(rownames(yfcst) == months_h[H]),,drop=FALSE]
      
        rownames(list_mean[[h]]) <- rownames(list_lower[[h]]) <- 
          rownames(list_upper[[h]]) <- months_h
      
        # see previous comments (first example)
        levels_h <- yfcst[months_h, target, drop=FALSE]
    
        ioae_h <- cbind(list_lower[[h]][,j], list_mean[[h]][,j],
                          list_upper[[h]][,j])
    
        colnames(ioae_h) <- c("lower", "mean", "upper")
        rownames(ioae_h) <- months_h
    
        
        Ty_h <- nrow(yfcst_h) - H
        yfcst_h[(Ty_h+1):(Ty_h+H), ] <- NA
        
        # see previous comments (first example)
        nowcast_h <- nowcasting_trans(ioae_h,
                            yfcst_h, target, trans, H)
      
        mae_levels_matrix[h, ] <- c(levels_h[j,], nowcast_h$LEVELS[j,])
      }
      mae_levels_list[[j]] <- mae_levels_matrix
    }
    
    # see previous comments (first example)
    MAE_levels_h <- colMeans(sapply(1:H, 
                                function(x) abs(mae_levels_list[[x]][,
                                "observed"] - mae_levels_list[[x]][,"mean"])))
    
    MedAE_levels_h <- apply(sapply(1:H, 
                                function(x) abs(mae_levels_list[[x]][,
                                "observed"] - mae_levels_list[[x]][,"mean"])),
                                2, median)
    
    # see previous comments (first example)
    mae_av_list <- mae_list
    MAE_av_h <- MAE_h
    MedAE_av_h <- MedAE
    
  }
  if(trans == "LEVELS"){
    # monthly variations
    mae_mv_list <- list()
    for(j in 1 : H){ 
      # see previous comments (first example)  
      mae_mv_matrix <- matrix(NA, Ht, 4)
      colnames(mae_mv_matrix) <- c("observed", "lower", "mean", "upper")
      rownames(mae_mv_matrix) <- rownames(xt)[out_sample[,j]]
      
      for(h in 1 : Ht){
        # see previous comments (first example)
        months_h <- rownames(x)[out_sample[h,]]
      
        yfcst_h <- yfcst[1:which(rownames(yfcst) == months_h[H]),,drop=FALSE]
      
        rownames(list_mean[[h]]) <- rownames(list_lower[[h]]) <- 
          rownames(list_upper[[h]]) <- months_h
      
        # see previous comments (first example)
        levels_h <- yfcst[months_h, target, drop=FALSE]
        mv_h <- cpm(yfcst_h[,target,drop=FALSE], 1)[months_h,,drop=FALSE]
    
      
        ioae_h <- cbind(list_lower[[h]][,j], list_mean[[h]][,j],
                          list_upper[[h]][,j])
    
        colnames(ioae_h) <- c("lower", "mean", "upper")
        rownames(ioae_h) <- months_h
    
        
        Ty_h <- nrow(yfcst_h) - H
        yfcst_h[(Ty_h+1):(Ty_h+H), ] <- NA
        
        # see previous comments (first example)
        nowcast_h <- nowcasting_trans(ioae_h,
                            yfcst_h, target, trans, H)
      
        mae_mv_matrix[h, ] <- c(mv_h[j,], nowcast_h$MV[j,])
      }
      mae_mv_list[[j]] <- mae_mv_matrix
    }
    
   # see previous comments (first example) 
    MAE_mv_h <- colMeans(sapply(1:H, 
                                function(x) abs(mae_mv_list[[x]][,"observed"] -
                                  mae_mv_list[[x]][,"mean"])))
    MedAE_mv_h <- apply(sapply(1:H, 
                                function(x) abs(mae_mv_list[[x]][,"observed"] -
                                  mae_mv_list[[x]][,"mean"])), 2, median)
    
    # annual variations
    mae_av_list <- list()
    for(j in 1 : H){ 
      # see previous comments (first example)  
      mae_av_matrix <- matrix(NA, Ht, 4)
      colnames(mae_av_matrix) <- c("observed", "lower", "mean", "upper")
      rownames(mae_av_matrix) <- rownames(xt)[out_sample[,j]]
      
      for(h in 1 : Ht){
       # see previous comments (first example) 
        months_h <- rownames(x)[out_sample[h,]]
      
        yfcst_h <- yfcst[1:which(rownames(yfcst) == months_h[H]),,drop=FALSE]
      
        rownames(list_mean[[h]]) <- rownames(list_lower[[h]]) <- 
          rownames(list_upper[[h]]) <- months_h
      
        # see previous comments (first example)
        levels_h <- yfcst[months_h, target, drop=FALSE]
        av_h <- cpm(yfcst_h[,target,drop=FALSE], 12)[months_h,,drop=FALSE]
    
      
        ioae_h <- cbind(list_lower[[h]][,j], list_mean[[h]][,j],
                          list_upper[[h]][,j])
    
        colnames(ioae_h) <- c("lower", "mean", "upper")
        rownames(ioae_h) <- months_h
    
        
        Ty_h <- nrow(yfcst_h) - H
        yfcst_h[(Ty_h+1):(Ty_h+H), ] <- NA
        
        # see previous comments (first example)
        nowcast_h <- nowcasting_trans(ioae_h,
                            yfcst_h, target, trans, H)
      
        mae_av_matrix[h, ] <- c(av_h[j,], nowcast_h$AV[j,])
      }
      mae_av_list[[j]] <- mae_av_matrix
    }
    
    # see previous comments (first example)
    MAE_av_h <- colMeans(sapply(1:H, 
                                function(x) abs(mae_av_list[[x]][,"observed"] -
                                  mae_av_list[[x]][,"mean"])))
    MedAE_av_h <- apply(sapply(1:H, 
                                function(x) abs(mae_av_list[[x]][,"observed"] -
                                  mae_av_list[[x]][,"mean"])), 2, median)
    
    # see previous comments (first example)
    mae_levels_list <- mae_list
    MAE_levels_h <- MAE_h
    MedAE_levels_h <- MedAE
  }
  # final result
  result <- list(mae_mv_list, MAE_mv_h, MedAE_mv_h,
                mae_av_list, MAE_av_h, MedAE_av_h,
                mae_levels_list, MAE_levels_h, MedAE_levels_h)
  
  names(result) <- c("mae_mv_list", "MAE_mv_h", "MedAE_mv_h",
                    "mae_av_list", "MAE_av_h", "MedAE_av_h",
                    "mae_levels_list", "MAE_levels_h", "MedAE_levels_h")
    
  return(result)
}

# @ pooled test
# e is the idiosyncratic errors
# type is the specification in the unit root test
# lag_max is the maximum number of lags to testing
pooled.test <- function(e, type = "const", lag_max = 7){
  N <- ncol(e)
  
  bics <- matrix(0, lag_max + 1, N)
  for(i in 1 : N){
    bics[,i] <- sapply(0:lag_max, function(x) adf(e[,i], type = type, 
                                                  k = x)$bic)
  }
  
  lags <- apply(bics, 2, function(x) which(x == min(x)) - 1)
  
  adf_ehat <- sapply(1:N, function(x) adf(e[,x], type = type,
                                          k = lags[x])$p.value)
  
  P <- (-2*sum(log(adf_ehat)) - 2*N)/sqrt(4*N)
  
  p.value <- dnorm(P)
  
  stats <- round(c(P, p.value), 4)
  names(stats) <- c("P", "p.value")
  
  result <- list(adf_ehat, stats)
  names(result) <- c("ADF", "stats")
  
  return(result)
}

