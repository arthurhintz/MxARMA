#' MxARMA fitting without regressors
#'
#' Used for fitting the Maxwell auto-regressive moving averages model without regressors. The model is given by
#'
#' \deqn{\eta_t = \alpha + \sum_{j\in ar}\phi_j \log(y_{t-j}) +\sum_{j\in ma}\theta_jr_{t-j}}
#' Where \itemize{
#' \item{\eqn{y} are the variables}
#' \item{\eqn{\alpha} is the intercept}
#' \item{\eqn{ar} are the indices for the auto-regression}
#' \item{\eqn{ma} are the indices for the moving-averages}
#' \item{\eqn{\phi} are the auto-regression coefficients}
#' \item{\eqn{\theta} are the moving-averages coefficients}
#' \item{\eqn{r} are the errors}}
#'
#'
#' @param y The vector of variables to fit
#' @param ar are the indices for the auto-regression.
#' @param ma are the indices for the moving-averages
#' @param h1 Number of predicted observations.
#' @param resid Residual type for the model (1 for quantile residuals, 2 for standardized residuals, 3 for alternative residuals).
#'
#'
#'
#'@examples
#'
#'y <- MxARMA::mxarma.sim(100, alpha=0.6, phi= c(0.6,0.2), theta = -0.4)
#'mxarma.fit(y, ar=c(1,2), ma=c(1), resid = 1, h1 = 12)
#'
#'
#'@importFrom stats make.link
#'@importFrom stats optim
#'@importFrom stats qnorm
#'@importFrom stats pnorm
#'
#'@export
mxarma.fit <- function(y, ar = NA, ma = NA, resid = 1, h1 = 0){

  if (!is.numeric(y) && min(y) < 0) {
    stop("y must be a numeric and positive vector")
  }

  isar <- !any(is.na(ar))
  isma <- !any(is.na(ma))

  link <- stats::make.link("log")
  linkfun <- link$linkfun
  linkinv <- link$linkinv

  p <- max(ar)
  q <- max(ma)
  n <- length(y)
  m <- max(p, q, na.rm=T)
  p1 <- length(ar)
  q1 <- length(ma)

  ylog = linkfun(y)

  y_prev <- c(rep(NA,(n+h1)))

  if (isar == T ){
  P <- matrix(rep(NA,(n-m)*p1),ncol=p1)

  for(i in 1:(n-m)){
    P[i,] <- ylog[i+m-ar]}

  Z <- cbind(rep(1,(n-m)), P)

  } else {
    Z <- as.matrix(rep(1,(n-m)))
  }

  # Sem regressores
  {x <- as.matrix(Z)
  Y <- y[(m+1):n]
  Ynew = linkfun(Y)
  ajuste = lm.fit(x, Ynew)
  mqo = c(ajuste$coef)}

#_____________________AR___________________________________
  if(isar == T && isma == F){

    loglik <- function(z){
      alpha <- z[1]
      phi = z[2:(p1+1)]

      eta <- rep(NA,n)

      for(i in (m+1):n){
        eta[i] <- alpha + (phi %*% ylog[i-ar]) }

      mu <- linkinv(eta[(m+1):n])
      y1 = y[(m+1):n]

      ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y1^2) + ((-4*y1^2)/(pi*mu^2))
      sum(ll)
    }

    escore <- function(z){
      alpha <- z[1]
      phi <- z[2:(p1+1)]

      eta<- rep(NA,n)

      for(i in (m+1):n){
        eta[i]<- alpha + (phi%*%ylog[i-ar]) }

      mu <- linkinv(eta[(m+1):n])
      y1 = y[(m+1):n]

      mT <- diag(mu)

      c <- ((8*y1^2)/(mu^3*pi))-(3/mu)

      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1,nrow=n)

      for(i in (m+1):n){
        deta.dalpha[i]<- 1
        deta.dphi[i,] <- P[(i-m),]
      }

      a <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]

      Ualpha <- t(a) %*% mT %*% c
      Uphi <-   t(rP) %*% mT %*% c

      return(c(Ualpha, Uphi))
  }

    opt <- stats::optim(mqo, fn = loglik, escore,
                 method = "BFGS",
                 control = list(fnscale = -1))


   names_par <- c("alpha", paste0("phi", ar))

   if (opt$conv != 0)
    warning("FUNCTION DID NOT CONVERGE!")

  z <- c()
  z$conv <- opt$conv
  coef <- opt$par
  names(coef)<-names_par
  z$coeff <- coef

  alpha <-coef[1]
  phi <- coef[2:(p1+1)]

  z$alpha <- alpha
  z$phi <- phi

  etahat<-rep(NA,n)
  errorhat<-rep(0,n)

  for(i in (m+1):n){
    etahat[i]<- alpha + (phi %*% ylog[i-ar])
    errorhat[i]<-ylog[i]-etahat[i]
  }
  muhat <- linkinv(etahat[(m+1):n])
  y1 = y[(m+1):n]

  z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
  z$etahat <- etahat
  z$errorhat <- errorhat


  deta.dalpha <- rep(0,n)
  deta.dphi <- matrix(0, ncol=p1,nrow=n)

  for(i in (m+1):n){
    deta.dalpha[i]<- 1
    deta.dphi[i,] <- P[(i-m),]
  }

  a <- deta.dalpha[(m+1):n]
  rP <- deta.dphi[(m+1):n,]

  W <- diag(((6)/(muhat^2))*(muhat^2))

  Kaa <- t(a) %*% W %*% a
  Kap <- t(a) %*% W %*% rP
  Kpa <- t(Kap)
  Kpp <- t(rP) %*% W %*% rP

  K <- rbind(
    cbind(Kaa,Kap),
    cbind(Kpa,Kpp)
  )

  z$K <- K

  #### Forecast

  ynew_prev <- c(ylog,rep(NA,h1))
  y_prev[1:n] <- z$fitted

  for(i in 1:h1){
    ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar])
    y_prev[n+i] <- linkinv(ynew_prev[n+i])
  }


}

#_________________MA_________________________________
  if(isar == F && isma == T){

    reg <- c(rep(0, q1 + 1))

    loglik <- function(z){
      alpha <- z[1]
      theta = z[2:(q1+1)]

      eta<-rep(NA,n)
      error <- rep(0,n)

      for(i in (m+1):n){
        eta[i] <- alpha + (theta %*% error[(i - ma)])
        error[i] <- ylog[i]-eta[i]}

      mu <- linkinv(eta[(m+1):n])
      y1 = y[(m+1):n]

      ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y1^2) + ((-4*y1^2)/(pi*mu^2))
      sum(ll)
    }

    escore <- function(z){
      alpha <- z[1]
      theta <- z[2:(q1+1)]

      eta <- rep(NA,n)
      error <- rep(0,n)

      for(i in (m+1):n){
        eta[i] <- alpha + (theta %*% error[(i - ma)])
        error[i] <- ylog[i]-eta[i]}

      mu <- linkinv(eta[(m+1):n])
      y1 = y[(m+1):n]

      mT <- diag(mu)

      c <- ((8*y1^2)/(mu^3*pi))-(3/mu)

      deta.dalpha <- rep(0,n)
      deta.dtheta<- matrix(0, ncol=q1,nrow=n)

      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)

      for(i in 1:(n-m)){
        R[i,] <- error[i+m-ma]
      }

      for(i in (m+1):n){
          deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
          deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
        }

      a <- deta.dalpha[(m+1):n]
      rR <- deta.dtheta[(m+1):n,]

      Ualpha <- t(a) %*% mT %*% c
      Utheta <-   t(rR) %*% mT %*% c

      return(c(Ualpha, Utheta))
    }

    opt <- stats::optim(reg, fn = loglik, escore,
                  method = "BFGS",
                  control = list(fnscale = -1))

    names_par <- c("alpha", paste0("theta", ma))

    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")

    z <- c()
    z$conv <- opt$conv
    coef <- opt$par
    names(coef)<-names_par
    z$coeff <- coef

    alpha <-coef[1]
    theta <- coef[2:(q1+1)]

    z$alpha <- alpha
    z$theta <- theta

    etahat<-rep(NA,n)
    errorhat<-rep(0,n)

    for(i in (m+1):n){
      etahat[i]<- alpha + (theta %*% errorhat[(i - ma)])
      errorhat[i]<- ylog[i]-etahat[i]}

    muhat <- linkinv(etahat[(m+1):n])

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat

    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m)){
        R[i,] <- errorhat[i+m-ma]}

    deta.dalpha <- rep(0,n)
    deta.dtheta <- matrix(0, ncol=q1,nrow=n)

    for(i in (m+1):n){
      deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
      deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
    }

    a <- deta.dalpha[(m+1):n]
    rR <- deta.dtheta[(m+1):n,]

    W <- diag(((6)/(muhat^2))*(muhat^2))

    Kaa <- t(a) %*% W %*% a
    Kat <- t(a) %*% W %*% rR
    Kta <- t(Kat)
    Ktt <- t(rR) %*% W %*% rR

    K <- rbind(
      cbind(Kaa,Kat),
      cbind(Kta,Ktt)
    )

    z$K <- K

    #### Forecasting
    ynew_prev <- c(ylog,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }

  }

#_________________________ARMA______________________
  if(isar == T && isma == T){

    reg <- c(mqo, rep(0,q1))

    loglik <- function(z){
      alpha <- z[1]
      phi = z[2:(p1+1)]
      theta = z[(p1+2):(q1+p1+1)]

      eta<-rep(NA,n)
      error<-rep(0,n)

      for(i in (m+1):n){
        eta[i] <- alpha + as.numeric(phi %*% ylog[i-ar]) + as.numeric(theta %*% error[i-ma])
        error[i] <- ylog[i]-eta[i]
        }

      mu <- linkinv(eta[(m+1):n])
      y1 = y[(m+1):n]

      ll <- log(32) - 2*log(pi) - 3*log(mu) + log(y1^2) + ((-4*y1^2)/(pi*mu^2))
      sum(ll)
    }

    escore <- function(z){
      alpha <- z[1]
      phi <- z[2:(p1+1)]
      theta <- z[(p1+2):(p1+q1+1)]

      error <- rep(0,n)
      eta <- rep(NA,n)

      for(i in (m+1):n){
        eta[i] <- alpha + as.numeric(phi %*% ylog[i - ar]) + as.numeric(theta %*% error[i - ma])
        error[i] <- ylog[i] - eta[i]
      }

      mu <- linkinv(eta[(m+1):n])
      y1 = y[(m+1):n]

      mT <- diag(mu)

      c <- ((8*y1^2)/(mu^3*pi))-(3/mu)

      R <- matrix(rep(NA,(n-m)*q1),ncol=q1)

      for(i in 1:(n-m)){
        R[i,] <- error[i+m-ma]
      }

      deta.dalpha <- rep(0,n)
      deta.dphi <- matrix(0, ncol=p1, nrow=n)
      deta.dtheta<- matrix(0, ncol=q1, nrow=n)

      for(i in (m+1):n){
        deta.dalpha[i]<- 1 - theta%*%deta.dalpha[i-ma]
        deta.dphi[i,]<- P[(i-m),] - theta%*%deta.dphi[i-ma,]
        deta.dtheta[i,]<- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
      }

      a <- deta.dalpha[(m+1):n]
      rP <- deta.dphi[(m+1):n,]
      rR <- deta.dtheta[(m+1):n,]

      Ualpha <- t(a) %*% mT %*% c
      Uphi <-   t(rP) %*% mT %*% c
      Utheta <- t(rR) %*% mT %*% c

      rval <- c(Ualpha,Uphi,Utheta)
    }


    opt <- stats::optim(reg,  loglik, escore,
                        method = "BFGS",
                        control = list(fnscale = -1))


    names_par <- c("alpha", paste0("phi", ar),  paste0("theta", ma))

    if (opt$conv != 0)
      warning("FUNCTION DID NOT CONVERGE!")

    z <- c()
    z$conv <- opt$conv
    coef <- opt$par
    names(coef) <- names_par
    z$coeff <- coef

    alpha <-coef[1]
    phi <- coef[2:(p1+1)]
    theta <- coef[(p1+2):(p1+q1+1)]

    z$alpha <- alpha
    z$phi <- phi
    z$theta <- theta

    errorhat<-rep(0,n)
    etahat<-rep(NA,n)

    for(i in (m+1):n){
      etahat[i] <- alpha + (phi %*% ylog[i-ar]) + (theta%*%errorhat[i-ma])
      errorhat[i] <- ylog[i]-etahat[i]}

    muhat <- linkinv(etahat[(m+1):n])
    y1 = y[(m+1):n]

    z$fitted <- ts(c(rep(NA,m),muhat),start=start(y),frequency=frequency(y))
    z$etahat <- etahat
    z$errorhat <- errorhat

    R <- matrix(rep(NA,(n-m)*q1),ncol=q1)
    for(i in 1:(n-m)){
      R[i,] <- errorhat[i+m-ma]
    }

    deta.dalpha <- rep(0,n)
    deta.dphi <- matrix(0, ncol=p1,nrow=n)
    deta.dtheta<- matrix(0, ncol=q1,nrow=n)

    for(i in (m+1):n){
      deta.dalpha[i] <- 1 - theta%*%deta.dalpha[i-ma]
      deta.dphi[i,] <- P[(i-m),] - theta%*%deta.dphi[i-ma,]
      deta.dtheta[i,] <- R[(i-m),] - theta%*%deta.dtheta[i-ma,]
    }

    a <- deta.dalpha[(m+1):n]
    rP <- deta.dphi[(m+1):n,]
    rR <- deta.dtheta[(m+1):n,]

    W <- diag(((6)/(muhat^2))*(muhat^2))

    Kaa <- t(a) %*% W %*% a
    Kap <- t(a) %*% W %*% rP
    Kpa <- t(Kap)
    Kat <- t(a) %*% W %*% rR
    Kta <- t(Kat)
    Kpp <- t(rP) %*% W %*% rP
    Kpt <- t(rP) %*% W %*% rR
    Ktp <- t(Kpt)
    Ktt <- t(rR) %*% W %*% rR

    K <- rbind(
      cbind(Kaa,Kap,Kat),
      cbind(Kpa,Kpp,Kpt),
      cbind(Kta,Ktp,Ktt))

    z$K <- K

    #### Forecasting
    ynew_prev <- c(ylog,rep(NA,h1))
    y_prev[1:n] <- z$fitted

    for(i in 1:h1){
      ynew_prev[n+i] <- alpha + (phi%*%ynew_prev[n+i-ar]) + (theta%*%errorhat[n+i-ma])
      y_prev[n+i] <- linkinv(ynew_prev[n+i])
      errorhat[n+i] <- 0
    }
  }

  ## Forecast

  z$y_prev = y_prev
  z$forecast <- y_prev[(n+1):(n+h1)]
  ############### Residuos

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

  ######################## Teste de hipótese

  vcov <- solve(K)
  z$vcov <- vcov

  stderror <- sqrt(diag(vcov))
  z$stderror <- stderror

  z$zstat <- abs(z$coeff/stderror)
  z$pvalues <- 2*(1 - stats::pnorm(z$zstat))

  z$loglik <- opt$value*(n/(n-m))
  z$counts <- as.numeric(opt$counts[1])


  z$aic <- -2*z$loglik+2*(p1+q1+2)
  z$bic <- -2*z$loglik+log(n)*(p1+q1+2)


  model_presentation <- cbind(round(z$coeff,4),round(z$stderror,4),round(z$zstat,4),round(z$pvalues,4))
  colnames(model_presentation)<-c("Estimate","Std. Error","z value","Pr(>|z|)")

  z$model <- model_presentation

  return(z)

}
