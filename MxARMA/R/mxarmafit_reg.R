#' MxARMA fitting with regressors
#'
#'Used for fitting the Maxwell auto-regressive moving averages model with regressors. The model is given by
#' \deqn{\eta_t = \beta_{0} + \boldsymbol{x}_t^\top\boldsymbol{\beta}+
#' \sum_{j\in ar}\phi_j [\log(y_{t-j})-\boldsymbol{x}_{t-j}^\top\boldsymbol{\beta}]
#' +\sum_{j\in ma}\theta_jr_{t-j}}
#'
#' #' Where \itemize{
#' \item{\eqn{y} are the variables}
#' \item{\eqn{\alpha} is the intercept}
#' \item{\eqn{\boldsymbol{x}} are the covariables}
#' \item{\eqn{\boldsymbol{\beta}} are the regression coefficients}
#' \item{\eqn{ar} are the indices for the auto-regression}
#' \item{\eqn{ma} are the indices for the moving-averages}
#' \item{\eqn{\phi} are the auto-regression coefficients}
#' \item{\eqn{\theta} are the moving-averages coefficients}
#' \item{\eqn{r} are the errors}}
#'
#'
#'
#' @param y The vector of variables to fit
#' @param cvar A matrix of covariates or regressors used in the model. It should have the appropriate dimensions for the computations.
#' @param ar Specified in the description.
#' @param ma Specified in the description.
#' @param h1 Number of predicted observations.
#' @param resid Residual type for the model (1 for quantile residuals, 2 for standardized residuals, 3 for alternative residuals).
#' @param X_hat covariate prediction vector
#'
#' @details The function fitted values from the MxARMA model, including the possibility of fitting values for an autoregressive with exogenous variables (ARX), moving average (MAX), or combined autoregressive moving average with exogenous variables (ARMAX) model.
#'
#' @examples
#'
#'  X = stats::runif(1004, 0, 1)
#'  y <- mxarmareg.sim(n = 1000, alpha = 0.6, phi = c(0.6, 0.1), theta = 0.3, beta = 0.7, X = X)
#'  X <- X[5:1004]
#'  xhat <- matrix(runif(12, 0, 1), ncol = 1)
#'  mxarmareg.fit(y = y, cvar = X, ar = c(1,2), ma = c(1), resid = 1, h1 = 12, X_hat = xhat)
#'
#'
#'@importFrom stats make.link
#'@importFrom stats optim
#'@importFrom stats qnorm
#'@importFrom stats pnorm
#'
#'
#' @export
mxarmareg.fit <- function(y, cvar, ar = NA, ma = NA, resid = 1, h1 = 0, X_hat = 0){

  if (!is.numeric(y) && min(y) < 0) {
    stop("y deve ser um vetor numérico e positivo")
  }

  isar <- !any(is.na(ar))
  isma <- !any(is.na(ma))

  y <- as.vector(y)
  X <- as.matrix(cvar)

  link <- make.link("log")
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  mu.eta <- link$mu.eta


  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p,q,na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)
  k1 <- ncol(X)

  ylog = linkfun(y)

  y_prev <- c(rep(NA,(n+h1)))

  if (isar == T ){
    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)

    for(i in 1:(n-m)){
      P[i,] <- ylog[i+m-ar]
    }

    Z <- cbind(rep(1,(n-m)),P)

  } else {
    Z <- as.matrix(rep(1,(n-m)))
  }

# com regressores
  {X_hat<-as.matrix(X_hat)
    x <- cbind(as.matrix(Z),X[(m+1):n,])
    Y <- y[(m+1):n]
    Ynew = linkfun(Y)
    ajuste = lm.fit(x, Ynew)
    mqo = c(ajuste$coef)}

#_____________ARX________________________
  if(isar == T && isma == F){
    print("ARX")

    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], beta1)

    loglik <- function(z){
      alpha <- z[1]
      phi <-  z[2:(p1+1)]
      beta <- z[(p1+2):length(z)]

      eta <- rep(NA,n)

      for(i in (m+1):n){
        eta[i] <-  alpha + X[i,] %*% as.matrix(beta) + (phi%*%(ylog[i-ar]-(X[i-ar,]%*%as.matrix(beta) )))
      }

      mu <- linkinv(eta[(m+1):n])

      y1 <-y[(m+1):n]

      ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y1^2) + ((-4*y1^2)/(pi*mu^2))
      sum(ll)
    }

    escore <- function(z){

      alpha <- z[1]
      phi <- z[2:(p1+1)]
      beta <- z[(p1+2):length(z)]

      eta<- rep(NA,n)

      for(i in (m+1):n){
        eta[i] <-  alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ylog[i-ar]- (X[i-ar,] %*%as.matrix(beta) )))
      }

      mu <- linkinv(eta[(m+1):n])
      y1 <-y[(m+1):n]

      P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
      for(i in 1:(n-m)){
        P[i,] <- ylog[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
      }

      M <- matrix(rep(NA,(n-m)*k1),ncol = k1)
      for(i in 1:(n-m)){
        M[i,] <- X[i + m,] - (phi %*% X[i + m - ar,])
      }

      mT <- diag(mu)

      c <- ((8*y1^2)/(mu^3*pi))-(3/mu)

      a <- matrix(rep(1,(n-m)),ncol=1)
      rP <- P
      rM <- M

      Ualpha <- t(a) %*% mT %*% c
      Uphi <-   t(rP) %*% mT %*% c
      Ubeta <-  t(rM) %*% mT %*% c

      rval <- c(Ualpha,Uphi,Ubeta)
    }

    opt <- optim(reg, loglik, escore,
                 method = "BFGS",
                 control = list(fnscale = -1))

    names_par <- c("alpha", paste0("phi", ar), paste0("beta", 1:k1))

    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")

    z <- c()
    z$conv <- opt$conv
    coef <- opt$par
    names(coef)<-names_par
    z$coeff <- coef

    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    beta <- coef[(p1+2):length(coef)]

    z$alpha <- alpha
    z$phi <- phi
    z$beta <- beta

    etahat<-rep(NA,n)

    for(i in (m+1):n){
      etahat[i] <-  alpha + X[i,] %*% as.matrix(beta) + (phi%*%(ylog[i-ar]-X[i-ar,]%*%as.matrix(beta) ))
    }

    muhat <- linkinv(etahat[(m+1):n])
    y1 <- y[(m+1):n]

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat


    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    for(i in 1:(n-m)){
      P[i,] <- ylog[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
    }

    M <- matrix(rep(NA,(n-m)*k1),ncol = k1)
    for(i in 1:(n-m)){
      M[i,] <- X[i + m,] - (phi %*% X[i + m - ar,])
    }

    a <- matrix(rep(1,(n-m)),ncol=1)
    rP <- P
    rM <- M

    W <- diag(((6)/(muhat^2))*(muhat^2))

    Kaa <- t(a) %*% W %*% a
    Kab <- t(a) %*% W %*% rM
    Kba <- t(Kab)
    Kap <- t(a) %*% W %*% rP
    Kpa <- t(Kap)
    Kbb <- t(rM) %*% W %*% rM
    Kbp <- t(rM) %*% W %*% rP
    Kpb <- t(Kbp)
    Kpp <- t(rP) %*% W %*% rP

    K <- rbind(
      cbind(Kaa,Kab,Kap),
      cbind(Kba,Kbb,Kbp),
      cbind(Kpa,Kpb,Kpp)
    )

    z$K <- K

    #### Forecasting
    ynew_prev <- c(ylog,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    X_prev<- rbind(X,X_hat)

    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) ))
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
    }

  }

#_____________MAX_____________________
  if(isar == F && isma == T) {
    print("MAX")

    beta1<- mqo[(2):length(mqo)]
    reg <- c(mqo[1], rep(0,q1), beta1)

    loglik <- function(z){
      alpha <- z[1]
      theta = z[(2):(q1+1)]
      beta <- z[(q1+2):length(z)]

      error <- rep(0,n)
      eta <- rep(NA,n)

      for(i in (m+1):n){
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        error[i] <- ylog[i]-eta[i]

      }
      mu <- linkinv(eta[(m+1):n])
      y1 <-  y[(m+1):n]

      ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y1^2) + ((-4*y1^2)/(pi*mu^2))
      sum(ll)
    }

    escore <- function(z){
      alpha <- z[1]
      theta <- z[(2):(q1+1)]
      beta <- z[(q1+2):length(z)]

      error <- rep(0,n)
      eta <- rep(NA,n)

      for(i in (m+1):n){
        eta[i]<- alpha + X[i,]%*%as.matrix(beta) + (theta%*%error[i-ma])
        error[i] <- ylog[i]-eta[i]
      }

      mu <- linkinv(eta[(m+1):n])
      y1 <-y[(m+1):n]

      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m)){
        R[i,] <- error[i+m-ma]
      }

      M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
      for(i in 1:(n-m)){
        for(j in 1:k1){
          M[i,j] <- X[i+m,j]
        }
      }

      deta.dalpha <- rep(0,n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)

      for(i in (m+1):n){
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
      }


      a <- deta.dalpha[(m+1):n]
      rR <- deta.dtheta[(m+1):n,]
      rM <- deta.dbeta[(m+1):n,]

      mT <- diag(mu)

      c <- ((8*y1^2)/(mu^3*pi))-(3/mu)

      Ualpha <- t(a) %*% mT %*% c
      Utheta <- t(rR) %*% mT %*% c
      Ubeta <-  t(rM) %*% mT %*% c

      rval <- c(Ualpha,Utheta,Ubeta)
    }

    opt <- optim(reg, loglik, escore,
                 method = "BFGS",
                 control = list(fnscale = -1))#, maxit = maxit1, reltol = 1e-12))

    names_par <- c("alpha", paste0("theta", ma), paste0("beta", 1:k1))

    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")

    z <- c()
    z$conv <- opt$conv
    coef <- opt$par
    names(coef)<-names_par
    z$coeff <- coef

    alpha <-coef[1]
    theta <- coef[2:(q1+1)]
    beta <- coef[(q1+2):length(coef)]

    z$alpha <- alpha
    z$theta <- theta
    z$beta <- beta

    errorhat <-rep(0,n)
    etahat <-rep(NA,n)

    for(i in (m+1):n){
      etahat[i] <- alpha + X[i,]%*%as.matrix(beta) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ylog[i]-etahat[i]
    }

    muhat <- linkinv(etahat[(m+1):n])
    y1 <- y[(m+1):n]

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat

    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m)){
      R[i,] <- errorhat[i+m-ma]
    }

    M <- matrix(rep(NA,(n-m)*k1),ncol=k1)
    for(i in 1:(n-m)){
      for(j in 1:k1){
        M[i,j] <- X[i+m,j]
      }
    }

    deta.dalpha <- rep(0,n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)

    for(i in (m+1):n){
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
    }

    a <- deta.dalpha[(m+1):n]
    rR <- deta.dtheta[(m+1):n,]
    rM <- deta.dbeta[(m+1):n,]

    W <- diag(((6)/(muhat^2))*(muhat^2))

    Kaa <- t(a) %*% W %*% a
    Kab <- t(a) %*% W %*% rM
    Kba <- t(Kab)
    Kat <- t(a) %*% W %*% rR
    Kta <- t(Kat)
    Kbb <- t(rM) %*% W %*% rM
    Kbt <- t(rM) %*% W %*% rR
    Ktb <- t(Kbt)
    Ktt <- t(rR) %*% W %*% rR

    K <- rbind(
      cbind(Kaa,Kab,Kat),
      cbind(Kba,Kbb,Kbt),
      cbind(Kta,Ktb,Ktt)
    )

    z$K <- K

# Forecasting
    ynew_prev <- c(ylog,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    X_prev<- rbind(X,X_hat)

    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }
  }

#_________________ARMAX_______________________________
  if(isar == T && isma == T){
    print("ARMAX")

    beta1<- mqo[(p1+2):length(mqo)]
    reg <- c(mqo[1:(p1+1)], rep(0,q1), beta1)

    loglik <- function(z){

      alpha <- z[1]
      phi = z[2:(p1+1)]
      theta = z[(p1+2):(p1+q1+1)]
      beta <- z[(p1+q1+2):length(z)]

      error <- rep(0,n)
      eta <- rep(NA,n)

      for(i in (m+1):n){
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ylog[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
        error[i] <- ylog[i]-eta[i]
      }

      mu <- linkinv(eta[(m+1):n])
      y1 <-  y[(m+1):n]

      ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y1^2) + ((-4*y1^2)/(pi*mu^2))
      sum(ll)
    }


    escore <- function(z){
      alpha <- z[1]
      phi <- z[2:(p1+1)]
      theta <- z[(p1+2):(p1+q1+1)]
      beta <- z[(p1+q1+2):length(z)]

      error<- rep(0,n)
      eta<- rep(NA,n)

      for(i in (m+1):n){
        eta[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ylog[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%error[i-ma])
        error[i] <- ylog[i]-eta[i]
      }

      mu <- linkinv(eta[(m+1):n])
      y1 <-y[(m+1):n]

      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
      for(i in 1:(n-m)){
        R[i,] <- error[i+m-ma]
      }

      P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
      for(i in 1:(n-m)){
        P[i,] <- ylog[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
      }

      M <- matrix(rep(NA,(n-m)*k1),ncol = k1)
      for(i in 1:(n-m)){
        M[i,] <- X[i + m,] - (phi %*% X[i + m - ar,])
      }

      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)
      deta.dbeta<- matrix(0, ncol=k1,nrow=n)

      for(i in (m+1):n){
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
      }

      a <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]
      rM <- deta.dbeta[(m+1):n,]

      mT <- diag(mu)
      c <- ((8*y1^2)/(mu^3*pi))-(3/mu)

      Ualpha <- t(a) %*% mT %*% c
      Uphi <-   t(rP) %*% mT %*% c
      Utheta <- t(rR) %*% mT %*% c
      Ubeta <-  t(rM) %*% mT %*% c

      rval <- c(Ualpha,Uphi,Utheta,Ubeta)
    }

    opt <- optim(reg, loglik, escore,
                 method = "BFGS",
                 control = list(fnscale = -1))

    names_par <- c("alpha", paste0("phi", ar), paste0("theta", ma), paste0("beta", 1:k1))

    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")

    z <- c()
    z$conv <- opt$conv
    coef <- opt$par
    names(coef)<-names_par
    z$coeff <- coef

    alpha <-coef[1]
    phi <-  coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]
    beta <- coef[(p1+q1+2):length(coef)]

    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta
    z$beta <- beta

    errorhat <-rep(0,n)
    etahat <-rep(NA,n)

    for(i in (m+1):n){
      etahat[i] <- alpha + X[i,]%*%as.matrix(beta) + (phi%*%(ylog[i-ar]-X[i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ylog[i] - etahat[i]
    }

    muhat <- linkinv(etahat[(m+1):n])
    y1 <- y[(m+1):n]

    z$errorhat <- errorhat

    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m)){
      R[i,] <- errorhat[i+m-ma]
    }

    P <- matrix(rep(NA,(n-m)*p1),ncol=p1)
    for(i in 1:(n-m)){
      P[i,] <- ylog[i+m-ar] - X[i+m-ar,]%*%as.matrix(beta)
    }

    M <- matrix(rep(NA,(n-m)*k1),ncol = k1)
    for(i in 1:(n-m)){
      M[i,] <- X[i + m,] - (phi %*% X[i + m - ar,])
    }

    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)
    deta.dbeta<- matrix(0, ncol=k1,nrow=n)

    for(i in (m+1):n)
    {
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      deta.dbeta[i,]<- M[(i-m),] - theta%*%deta.dbeta[i-ma,]
    }

    a <- deta.dalpha[(m+1):n]
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]
    rM <- deta.dbeta[(m+1):n,]


    W<- diag(((6)/(muhat^2))*(muhat^2))

    Kaa <- t(a) %*% W %*% a
    Kab <- t(a) %*% W %*% rM
    Kba <- t(Kab)
    Kap <- t(a) %*% W %*% rP
    Kpa <- t(Kap)
    Kat <- t(a) %*% W %*% rR
    Kta <- t(Kat)

    Kbb <- t(rM) %*% W %*% rM
    Kbp <- t(rM) %*% W %*% rP
    Kpb <- t(Kbp)
    Kbt <- t(rM) %*% W %*% rR
    Ktb <- t(Kbt)

    Kpp <- t(rP) %*% W %*% rP
    Kpt <- t(rP) %*% W %*% rR
    Ktp <- t(Kpt)

    Ktt <- t(rR) %*% W %*% rR

    K <- rbind(
      cbind(Kaa,Kap,Kat,Kab),
      cbind(Kpa,Kpp,Kpt,Kpb),
      cbind(Kta,Ktp,Ktt,Ktb),
      cbind(Kba,Kbp,Kbt,Kbb)
    )

    z$K <- K

    fited <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$fitted = fited
    z$etahat <- etahat

#Forecasting
    ynew_prev <- c(ylog,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    X_prev<- rbind(X,X_hat)

    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + X_prev[n+i,]%*%as.matrix(beta) + (phi%*%(ynew_prev[n+i-ar]-X_prev[n+i-ar,]%*%as.matrix(beta) )) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }

  }

#Forecast

  z$y_prev = y_prev
  z$forecast <- y_prev[(n+1):(n+h1)]


#Residuos

var_y <- as.vector((z$fitted^2)*(0.178))

resid1 <- as.vector(stats::qnorm(MxARMA::pmax(y,z$fitted))) # quantile residuals

resid2 <- as.vector((y-z$fitted)/sqrt(var_y)) # standardized

l_tilde <- log(MxARMA::dmax(y[(m+1):n], y[(m+1):n]))
l_hat <- log(MxARMA::dmax(y[(m+1):n], z$fitted[(m+1):n] ))

dt <- abs(l_tilde-l_hat)

resid3 <- as.vector(c(rep(NA,m),sign(y[(m+1):n]-z$fitted[(m+1):n])*sqrt(2*(dt))))

if(resid==1) residc <- resid1[(m+1):n]
if(resid==2) residc <- resid2[(m+1):n]
if(resid==3) residc <- resid3[(m+1):n]

z$residc <- residc

#Teste de hipótese

vcov <- solve(K)
z$vcov <- vcov

stderror <- sqrt(diag(vcov))
z$stderror <- stderror

z$zstat <- abs(z$coeff/stderror)
z$pvalues <- 2*(1 - stats::pnorm(z$zstat))


z$loglik <- opt$value*(n/(n - m))
z$counts <- as.numeric(opt$counts[1])

z$aic <- -2*z$loglik+2*(p1+q1+2+length(beta))
z$bic <- -2*z$loglik+log(n)*(p1+q1+2+length(beta))


model_presentation <- cbind(round(z$coeff,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
colnames(model_presentation) <- c("Estimate","Std. Error","z value","Pr(>|z|)")

z$model <- model_presentation

return(z)

}
