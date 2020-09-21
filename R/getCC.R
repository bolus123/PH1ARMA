simInnov <- function(n, sigma = 1, XSim = 'norm', XPars = c(0, 1)) {
  
  if (XSim == "norm") {
    
    me <- XPars[1];
    std <- sqrt(XPars[2]);
    
    pars <- c(me, std)
    
    rgen <- rnorm
    
  } else if (XSim == "t") {
    
    me <- 0;
    std <- sqrt(XPars[1] / (XPars[1] - 2));
    
    pars <- c(XPars[1])
    
    rgen <- rt
    
  } else if (XSim == "gamma") {
    
    me <- XPars[1] * XPars[2];
    std <- sqrt(XPars[1] * XPars[2] ^ 2);
    
    pars <- c(XPars[1], 1 / XPars[2])
    
    rgen <- rgamma
    
    
  } else if (XSim == "beta") {
    
    me <- XPars[1] / (XPars[1] + XPars[2]);
    std <- sqrt(XPars[1] * XPars[2] / ((XPars[1] + XPars[2]) ^ 2) / (XPars[1] + XPars[2] + 1));
    
    pars <- c(XPars[1], XPars[2])
    
    rgen <- rbeta
    
  }
  
  if (XSim != 't') {
    
    out <- rgen(n, pars[1], pars[2])
    
  } else {
    
    out <- rgen(n, pars[1])
    
  }
  
  (out - me) / std * sigma
  
}

simARMAProcess <- function(n, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, sigma = 1, innovDist = 'norm', innovPars = c(0, 1), burnIn = 1000) {
  
  if (burnIn > 0) {
    
    arima.sim(list(order = order, ar = phiVec, ma = thetaVec),
              n = n + burnIn,
              innov = simInnov(n + burnIn, sigma = sigma, XSim = innovDist, XPars = innovPars))[-(1:burnIn)]
    
  } else {
    
    arima.sim(list(order = order, ar = phiVec, ma = thetaVec),
              n = n + burnIn,
              innov = simInnov(n + burnIn, sigma = sigma, XSim = innovDist, XPars = innovPars))
    
  }
  
  
}



simCoefDist <- function(n, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, method = 'CSS-ML',
                        nsim = 1000, burnIn = 1000, seed = 12345) {
  
  set.seed(seed)
  
  outAR <- matrix(NA, nrow = nsim, ncol = length(phiVec))
  outMA <- matrix(NA, nrow = nsim, ncol = length(thetaVec))
  outMean <- rep(NA, nsim)
  outGamma <- rep(NA, nsim)
  
  for (i in 1:nsim) {
    
      flg <- 1
    
      while (flg == 1) {
      
        sim <- simARMAProcess(n, order, phiVec, thetaVec, sigma = 1, innovDist = 'norm', innovPars = c(0, 1), burnIn = burnIn)
        model <- try(arima(sim, order = order, method = method), silent = TRUE)
        
        if (class(model) != 'try-error') {
          if (!is.null(phiVec)) {
            outAR[i, ] <- model$coef[1:order[1]]
          } else {
            outAR[i, ] <- NULL
          }
          
          if (!is.null(thetaVec)) {
            outMA[i, ] <- model$coef[(order[1] + 1):(order[1] + order[3])]
          } else {
            outMA[i, ] <- NULL
          }
          
          outMean[i] <- mean(sim)
          outGamma[i] <- var(sim)
          
          flg <- 0
          
        } 
        
      }
    
    
  }
  
  list(phiVec = outAR, thetaVec = outMA)
  
}


fapPH1ARMAUnknown <- function(cc = 3, n = 30, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL,
                       nsim = 1000, burnIn = 1000) {
  
  out <- lapply(1:nsim, function(X){
    sim <- simARMAProcess(n, order, phiVec, thetaVec, sigma = 1,
                          innovDist = 'norm', innovPars = c(0, 1), burnIn = burnIn)
    sim <- (sim - mean(sim)) / sd(sim)
    sum(-cc <= sim & sim <= cc) != n
  })
  
  mean(unlist(out))
  
}

fapPH1ARMAKnown <- function(cc = 3, n = 30, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, gamma0 = 1,
                       nsim = 1000, burnIn = 1000) {
  
  #cat('phi:', phi, '\n')
  #cat('gamma0:', gamma0, '\n')
  out <- lapply(1:nsim, function(X){
    sim <- simARMAProcess(n, order, phiVec, thetaVec, sigma = 1,
                          innovDist = 'norm', innovPars = c(0, 1), burnIn = burnIn)
    sim <- sim / sqrt(gamma0)
    sum(-cc <= sim & sim <= cc) != n
  })
  
  mean(unlist(out))
  
}

getCCPH1ARMA <- function(FAP0 = 0.1, interval = c(1, 4), n = 30, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL,
                         nsim = 1000, burnIn = 1000, seed = 12345) {
  
  root.finding <- function(FAP0, cc, n, order, phiVec, thetaVec, nsim1, nsim2, burnIn, seed) {
    
    set.seed(seed)
    
    if (nsim1 > 0) {
      #cat('Unknown', '\n')
      FAPin <- lapply(1:nsim1, function(X) {
        
        phiVecTmp <- phiVec[X, ]
        if (length(phiVecTmp) == 0) phiVecTmp <- NULL
        
        thetaVecTmp <- thetaVec[X, ]
        if (length(thetaVecTmp) == 0) thetaVecTmp <- NULL
        
        fapPH1ARMAUnknown(cc, n, order, phiVecTmp, thetaVecTmp, nsim2, burnIn)
        
      })
      
      FAPin <- mean(unlist(FAPin))
      
    } else {
      #cat('Known', '\n')
	  if (order == c(1, 0, 0)) {
		gamma0 <- 1 / (1 - phiVec ^ 2)
	  } else if (order == c(0, 0, 1)) {
		gamma0 <- (1 + thetaVec ^ 2)
	  } else if (order == c(1, 0, 1)) {
		gamma0 <- (1 + 2 * phiVec * thetaVec + thetaVec ^ 2) / (1 - phiVec ^ 2)
	  } else {
		sim <- simARMAProcess(50000000, order, phiVec, thetaVec, sigma = 1,
                          innovDist = 'norm', innovPars = c(0, 1), burnIn = burnIn)
		gamma0 <- var(sim)
	  }
	 
      FAPin <- fapPH1ARMAKnown(cc, n, order, phiVec, thetaVec, gamma0,
                       nsim2, burnIn)
      
    }
    
    cat('FAPin:', FAPin, 'and cc:', cc, '\n')
    FAP0 - FAPin
    
  }
  
  if (is.matrix(phiVec) | is.matrix(thetaVec)) {
    nsim1 <- max(dim(phiVec)[1], dim(thetaVec)[1])
  } else {
    nsim1 <- 0
  }
  
  uniroot(root.finding, interval, FAP0 = FAP0, n = n, order = order, phiVec = phiVec, thetaVec = thetaVec,
          nsim1 = nsim1, nsim2 = nsim, burnIn = burnIn, seed = seed)$root
  
  
}

getCC <- function(FAP0 = 0.1, interval = c(1, 4), n = 30, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, method = 'CSS-ML',
                  double.sim = TRUE, nsimCoefs = 1000, nsimProcess = 1000, burnIn = 1000, seed = 12345) {
  
  
  if (double.sim == TRUE) {
    
    CoefDist <- simCoefDist(n, order, phiVec, thetaVec, method, nsim = nsimCoefs, burnIn = burnIn, seed = seed)
    
    out <- getCCPH1ARMA(FAP0, interval, n, order, phiVec = CoefDist$phiVec, thetaVec = CoefDist$thetaVec,
                        nsim = nsimProcess, burnIn = burnIn, seed = seed)
    
  } else {
    
    out <- getCCPH1ARMA(FAP0, interval, n, order, phiVec = phiVec, thetaVec = thetaVec,
                        nsim = nsimProcess, burnIn = burnIn, seed = seed)
    
  }
  
  out
  
}


PH1ARMA <- function(X, cc = NULL, FAP0 = 0.1, order = NULL, method = 'CSS-ML', plot.option = TRUE,
                    interval = c(1, 4), double.sim = TRUE, nsimCoefs = 1000, nsimProcess = 1000, burnIn = 1000, seed = 12345) {
  
  Y <- X
  
  n <- length(Y)
  
  if (is.null(cc)) {
    
    if (is.null(order)) {
      
      model <- auto.arima(Y, method = method)
      order <- rep(0, 3)
      order[1] <- length(model$model$phi)
      order[2] <- length(model$model$Delta)
      order[3] <- length(model$model$theta)
      
    } else {
      
      model <- arima(Y, order = order, method = method)
      
    }
    
    if (length(model$model$phi) > 0) {
      phiVec <- model$model$phi
    } else {
      phiVec <- NULL
    }
    
    if (length(model$model$theta) > 0) {
      thetaVec <- model$model$theta
    } else {
      thetaVec <- NULL
    }
    
    cc <- getCC(FAP0 = FAP0, interval = interval, n, order = order, phiVec = phiVec, thetaVec = thetaVec, method = method,
                double.sim = double.sim, nsimCoefs = nsimCoefs, nsimProcess = nsimProcess, burnIn = burnIn, seed = seed)
    
  }
  
  if (order[2] > 0) {
    
    Y <- diff(Y, differences = order[2])
    
  }
  
  mu <- mean(Y)
  gamma <- sd(Y)
  
  stdX <- (Y - mu) / gamma
  
  LCL <- -cc
  UCL <- cc
  
  if (plot.option == TRUE) {
    
    main.text <- paste('Phase I Individual Chart for FAP0 =', FAP0, 'with an ARMA model')
    
    plot(c(1, n), c(min(LCL, stdX), max(UCL, stdX)), xaxt = "n", xlab = 'Observation', ylab = 'Charting Statistic', type = 'n', main = main.text)
    
    axis(side = 1, at = 1:n)
    
    points(1:n, stdX, type = 'o')
    points(c(-1, n + 2), c(LCL, LCL), type = 'l', col = 'red')
    points(c(-1, n + 2), c(UCL, UCL), type = 'l', col = 'red')
    points(c(-1, n + 2), c(mu, mu), type = 'l', col = 'blue')
    text(round(n * 0.8), UCL, paste('UCL = ', round(UCL, 4)), pos = 1)
    text(round(n * 0.8), LCL, paste('LCL = ', round(LCL, 4)), pos = 3)
    
  }
  
  list(CL = mu, gamma = gamma, cc = cc, order = order, phiVec = phiVec, thetaVec = thetaVec, LCL = LCL, UCL = UCL, CS = stdX)
  
  
}
