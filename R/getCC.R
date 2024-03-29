##############################################
  #Matrix form
##############################################
InvertQ <- function(coef){
  
  out <- 1
  
  if (all(abs(coef) < 1)) {
    
    minmod <- min(Mod(polyroot(c(1, coef))))
    
    if (minmod <= 1) {
      
      out <- 0
      
    }
    
  } else {
    
    out <- 0
    
  }

  return(out)
  
}

parsMat <- function(n, parsVec, norder = 1) {
  Check <- InvertQ(parsVec)
  if (Check == 0) {
    NULL
  } else {
    Mat <- diag(n)
    for (i in 1:norder) {
      Mat <- Mat + Diag(rep(parsVec[i], n - i), k = -i)
    }
    Mat
  }
}



SigmaMat <- function(n, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, sigma2 = 1, burnIn = 50) {
  if (order[1] == 0) {
    phiMat <- diag(n + burnIn)
  } else {
    phiMat <- parsMat(n + burnIn, -phiVec, norder = order[1])
  }

  if (order[3] == 0) {
    thetaMat <- diag(n + burnIn)
  } else {
    thetaMat <- parsMat(n + burnIn, thetaVec, norder = order[3])
  }

  out <- solve(phiMat) %*% thetaMat %*% t(thetaMat) %*% t(solve(phiMat)) * sigma2

  gamma0 <- out[dim(out)[1], dim(out)[2]]

  if (burnIn > 0) {
    out <- out[-c(1:burnIn), -c(1:burnIn)]
  }

  list(SigmaMat = out, sqrtSigmaMat = sqrtm(out)$B, gamma0 = gamma0)
}

##############################################
  #process simulation
##############################################

# Simulate innovations
simInnov <- function(n, sigma2 = 1, XSim = 'norm', XPars = c(0, 1)) {

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

  (out - me) / std * sqrt(sigma2)

}



simARMAProcess <- function(n, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, sigma2 = 1, innovDist = 'norm', innovPars = c(0, 1), burnIn = 50,
                           simType = 'Matrix', SigMat = SigmaMat(n = n, order = order, phiVec = phiVec, thetaVec = thetaVec, sigma2 = sigma2, burnIn = burnIn)) {

  if (simType == 'Matrix') {

    as.vector(SigMat$sqrtSigmaMat %*% matrix(simInnov(n, sigma2 = 1, XSim = innovDist, XPars = innovPars), ncol = 1))

  } else if (simType == 'Recursive') {

    innov <- simInnov(n + burnIn + order[1] + order[3], sigma2 = sigma2, XSim = innovDist, XPars = innovPars)

    n.start <- order[1] + order[3]

    if (burnIn > 0) {

      arima.sim(list(order = order, ar = phiVec, ma = thetaVec),
                n = n + burnIn, innov = innov[-c(1:(n.start))], n.start = n.start, start.innov = innov[c(1:(n.start))])[(burnIn + 1):(n + burnIn)]

    } else {

      arima.sim(list(order = order, ar = phiVec, ma = thetaVec),
                n = n + burnIn, innov = innov[-c(1:(n.start))], n.start = n.start, start.innov = innov[c(1:(n.start))])

    }

  }


}


###simARMAProcessMissingValue <- function(n, mvVec = NULL, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, sigma2 = 1, innovDist = 'norm', innovPars = c(0, 1), burnIn = 50,
###                           simType = 'Matrix', SigMat = SigmaMat(n = n, order = order, phiVec = phiVec, thetaVec = thetaVec, sigma2 = sigma2, burnIn = burnIn)) {
###
###  if (simType == 'Matrix') {
###
###    as.vector(SigMat$sqrtSigmaMat %*% matrix(simInnov(n, sigma2 = 1, XSim = innovDist, XPars = innovPars), ncol = 1))
###
###  } else if (simType == 'Recursive') {
###
###    innov <- simInnov(n + burnIn + order[1] + order[3], sigma2 = sigma2, XSim = innovDist, XPars = innovPars)
###
###    n.start <- order[1] + order[3]
###
###    if (burnIn > 0) {
###
###      sim <- arima.sim(list(order = order, ar = phiVec, ma = thetaVec),
###                n = n + burnIn, innov = innov[-c(1:(n.start))], n.start = n.start, start.innov = innov[c(1:(n.start))])[(burnIn + 1):(n + burnIn)]
###
###    } else {
###
###      sim <- arima.sim(list(order = order, ar = phiVec, ma = thetaVec),
###                n = n + burnIn, innov = innov[-c(1:(n.start))], n.start = n.start, start.innov = innov[c(1:(n.start))])
###
###    }
###
###  }
###
###  if (max(mvVec) <= n) {
###
###    mv <- rep(0, n)
###    mv[mvVec] <- 1
###
###  }
###
###  list(sim = sim, mv = mv)
###
###}


simCoefDist <- function(n, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, method = 'Method 3',
                        nsim = 100, burnIn = 50, simType = 'Matrix',
                        SigMat = SigmaMat(n, order, phiVec, thetaVec, sigma2 = 1, burnIn = burnIn), seed = 12345) {

  set.seed(seed)

  outAR <- matrix(NA, nrow = nsim, ncol = order[1])
  outMA <- matrix(NA, nrow = nsim, ncol = order[3])
  outMean <- rep(NA, nsim)
  outGamma <- rep(NA, nsim)

  for (i in 1:nsim) {

      flg <- 1

      while (flg == 1) {

        sim <- simARMAProcess(n, order, phiVec, thetaVec, sigma2 = 1, innovDist = 'norm', innovPars = c(0, 1), burnIn = burnIn,
                              simType = simType, SigMat = SigMat)

        if (method == 'Method 1' | method == 'Method 3') {
          model <- try(arima(sim, order = order, method = 'CSS-ML'), silent = TRUE)
        } else if (method == 'Method 2') {
          model <- try(arima(sim, order = order, method = 'CSS'), silent = TRUE)
        }

        check1 <- 1
        check2 <- 1

        if (class(model) != 'try-error') {
          if (order[1] > 0) {
            outAR[i, ] <- model$coef[1:order[1]]
            check1 = InvertQ(outAR[i, ]) == 1
            #cat('AR:', outAR[i, ], '\n')
            #cat('check1:', check1, '\n')
          } else {
            outAR[i, ] <- rep(0, order[1])
          }

          if (order[3] > 0) {
            outMA[i, ] <- model$coef[(order[1] + 1):(order[1] + order[3])]
            check2 = InvertQ(outMA[i, ]) == 1
            #cat('MA:', outMA[i, ], '\n')
            #cat('check2:', check2, '\n')
          } else {
            outMA[i, ] <- rep(0, order[3])
          }


          if (check1 & check2) {

            outMean[i] <- mean(sim)
            outGamma[i] <- var(sim)

            flg <- 0
          }

        }

      }


  }

  list(phiVec = outAR, thetaVec = outMA)

}


simCoefDistMissingValue <- function(n, mvVec = NULL, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, method = 'Method 3',
                        nsim = 100, burnIn = 50, simType = 'Matrix',
                        SigMat = SigmaMat(n, order, phiVec, thetaVec, sigma2 = 1, burnIn = burnIn),
                        tol = 1e-2, seed = 12345) {

  set.seed(seed)

  outAR <- matrix(NA, nrow = nsim, ncol = order[1])
  outMA <- matrix(NA, nrow = nsim, ncol = order[3])
  outMean <- rep(NA, nsim)
  outGamma <- rep(NA, nsim)

  mv <- rep(0, n)
  mv[mvVec] <- 1

  for (i in 1:nsim) {

    flg <- 1

    while (flg == 1) {

      sim <- simARMAProcess(n, order, phiVec, thetaVec, sigma2 = 1, innovDist = 'norm', innovPars = c(0, 1), burnIn = burnIn,
                            simType = simType, SigMat = SigMat)

      oldLogLikelihood <- -Inf
      newLogLikelihood <- Inf

      while(abs(newLogLikelihood - oldLogLikelihood) > tol) {

        if (method == 'Method 1' | method == 'Method 3') {
          model <- try(arima(sim, order = order, method = 'CSS-ML'), silent = TRUE)
        } else if (method == 'Method 2') {
          model <- try(arima(sim, order = order, method = 'CSS'), silent = TRUE)
        }

        if (class(model) != 'try-error') {

          sim[mv] <- sim[mv] - model$residuals[mv]
          oldLogLikelihood <- newLogLikelihood
          newLogLikelihood <- model$loglik

        } else {

          newLogLikelihood <- 0
          oldLogLikelihood <- 0

        }

        #cat("new:", newLogLikelihood, 'and old:', oldLogLikelihood, '\n')
        #cat('model', model$coef, '\n')

      }

      check1 <- 1
      check2 <- 1

      if (class(model) != 'try-error') {
        if (order[1] > 0) {
          outAR[i, ] <- model$coef[1:order[1]]
          check1 = InvertQ(outAR[i, ]) == 1
          #cat('AR:', outAR[i, ], '\n')
          #cat('check1:', check1, '\n')
        } else {
          outAR[i, ] <- rep(0, order[1])
        }

        if (order[3] > 0) {
          outMA[i, ] <- model$coef[(order[1] + 1):(order[1] + order[3])]
          check2 = InvertQ(outMA[i, ]) == 1
          #cat('MA:', outMA[i, ], '\n')
          #cat('check2:', check2, '\n')
        } else {
          outMA[i, ] <- rep(0, order[3])
        }


        if (check1 & check2) {

          outMean[i] <- mean(sim)
          outGamma[i] <- var(sim)

          flg <- 0
        }

      }

    }


  }

  list(phiVec = outAR, thetaVec = outMA)

}



fapPH1ARMA <- function(cc = 3, n = 50, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, case = 'U', method = 'Method 3',
                       nsim = 1000, burnIn = 50, simType = 'Matrix') {

  if (simType == 'Matrix') {
      sigMat1 <- SigmaMat(n, order, phiVec, thetaVec, sigma2 = 1, burnIn = burnIn)
      gamma0 <- sigMat1$gamma0
  } else if (simType == 'Recursive') {
    if (case == 'K') {
      sigMat1 <- SigmaMat(100, order, phiVec, thetaVec, sigma2 = 1, burnIn = 50)
      gamma0 <- sigMat1$gamma0
    } else {
      sigMat1 <- NULL
    }
  }

  out <- lapply(1:nsim, function(X){
    flg <- 1
    while (flg == 1) {
      sim <- simARMAProcess(n, order, phiVec, thetaVec, sigma2 = 1,
                            innovDist = 'norm', innovPars = c(0, 1), burnIn = burnIn, simType = simType, SigMat = sigMat1)

      if (case == 'U') {
        if (method == 'Method 2' | method == 'Method 3') {

          sim <- (sim - mean(sim)) / sd(sim)
          flg <- 0
          out1 <- sum(-cc <= sim & sim <= cc) != n
          return(out1)

        } else if (method == 'Method 1') {

          m1 <- try(arima(sim, order = order, method = 'CSS-ML'), silent = TRUE)

          if (class(m1) != 'try-error') {

            if (!is.null(phiVec)) {
              phiVecNew <- m1$coef[1:order[1]]
            } else {
              phiVecNew <- NULL
            }

            if (!is.null(thetaVec)) {
              thetaVecNew <- m1$coef[(order[1] + 1):(order[1] + order[3])]
            } else {
              thetaVecNew <- NULL
            }

            Intercept <- m1$coef[order[1] + order[3] + 1]
            mu0 <- Intercept * (1 - sum(phiVecNew))
            sigma2 <- m1$sigma2
            gamma0 <- SigmaMat(100, order = order, phiVec = phiVec, thetaVec = thetaVec, sigma2 = sigma2, burnIn = 50)$gamma0

            sim <- (sim - mu0) / sqrt(gamma0)
            flg0 <- 0
            out1 <- sum(-cc <= sim & sim <= cc) != n
            return(out1)

          }
        }
  	  } else if (case == 'K') {

  		  sim <- (sim) / sqrt(gamma0)
  		  flg <- 0
  		  out1 <- sum(-cc <= sim & sim <= cc) != n
  		  return(out1)

  	  }


    }
  })

  mean(unlist(out))

}


fapPH1ARMAMissingValue <- function(cc = 3, n = 50, mvVec = NULL, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, case = 'U', method = 'Method 3',
                       nsim = 1000, burnIn = 50, simType = 'Matrix', tol = 1e-2) {

  if (simType == 'Matrix') {
    sigMat1 <- SigmaMat(n, order, phiVec, thetaVec, sigma2 = 1, burnIn = burnIn)
    gamma0 <- sigMat1$gamma0
  } else if (simType == 'Recursive') {
    if (case == 'K') {
      sigMat1 <- SigmaMat(100, order, phiVec, thetaVec, sigma2 = 1, burnIn = 50)
      gamma0 <- sigMat1$gamma0
    } else {
      sigMat1 <- NULL
    }
  }

  mv <- rep(0, n)
  mv[mvVec] <- 1

  out <- lapply(1:nsim, function(X){
    flg <- 1
    while (flg == 1) {
      sim <- simARMAProcess(n, order, phiVec, thetaVec, sigma2 = 1,
                            innovDist = 'norm', innovPars = c(0, 1), burnIn = burnIn, simType = simType, SigMat = sigMat1)


      if (case == 'U') {
        if (method == 'Method 2' | method == 'Method 3') {

          sim <- sim[-mvVec]

          sim <- (sim - mean(sim)) / sd(sim)
          flg <- 0
          out1 <- sum(-cc <= sim & sim <= cc) != n - length(mvVec)
          return(out1)

        } else if (method == 'Method 1') {

          oldLogLikelihood <- -Inf
          newLogLikelihood <- Inf

          while(abs(newLogLikelihood - oldLogLikelihood) > tol) {

            m1 <- try(arima(sim, order = order, method = 'CSS-ML'), silent = TRUE)

            if (class(m1) != 'try-error') {

              sim[mv] <- sim[mv] - m1$residuals[mv]
              oldLogLikelihood <- newLogLikelihood
              newLogLikelihood <- m1$loglik

            } else {

              newLogLikelihood <- 0
              oldLogLikelihood <- 0

            }

          }

          if (class(m1) != 'try-error') {

            if (!is.null(phiVec)) {
              phiVecNew <- m1$coef[1:order[1]]
            } else {
              phiVecNew <- NULL
            }

            if (!is.null(thetaVec)) {
              thetaVecNew <- m1$coef[(order[1] + 1):(order[1] + order[3])]
            } else {
              thetaVecNew <- NULL
            }

            Intercept <- m1$coef[order[1] + order[3] + 1]
            mu0 <- Intercept * (1 - sum(phiVecNew))
            sigma2 <- m1$sigma2
            gamma0 <- SigmaMat(100, order = order, phiVec = phiVec, thetaVec = thetaVec, sigma2 = sigma2, burnIn = 50)$gamma0

            sim <- (sim - mu0) / sqrt(gamma0)
            flg0 <- 0
            out1 <- sum(-cc <= sim & sim <= cc) != n
            return(out1)

          }
        }
      } else if (case == 'K') {

        sim <- sim[-mvVec]

        sim <- (sim) / sqrt(gamma0)
        flg <- 0
        out1 <- sum(-cc <= sim & sim <= cc) != n
        return(out1)

      }


    }
  })

  mean(unlist(out))

}


getCCPH1ARMASim <- function(FAP0 = 0.1, interval = c(1, 4), n = 50, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, case = 'U', method = 'Method 3',
                         nsim = 1000, burnIn = 50, simType = 'Matrix', seed = 12345) {

  root.finding <- function(FAP0, cc, n, order, phiVec, thetaVec, case, method, nsim1, nsim2, burnIn, simType, seed) {

    set.seed(seed)

    if (nsim1 > 0) {

      FAPin <- lapply(1:nsim1, function(X) {

        phiVecTmp <- phiVec[X, ]
        if (length(phiVecTmp) == 0) phiVecTmp <- rep(0, order[1])

        thetaVecTmp <- thetaVec[X, ]
        if (length(thetaVecTmp) == 0) thetaVecTmp <- rep(0, order[3])

        #cat('AR:', phiVecTmp, '\n')
        #cat('checkAR:', all(abs(outAR[i, ]) < 1), '\n')

        #cat('MA:', thetaVecTmp, '\n')
        #cat('checkMA:', all(abs(outMA[i, ]) < 1), '\n')

        fapPH1ARMA(cc = cc, n = n, order = order, phiVec = phiVecTmp, thetaVec = thetaVecTmp,
			    case = case, method = method, nsim = nsim2, burnIn = burnIn, simType = simType)

      })

      FAPin <- mean(unlist(FAPin))

    } else {


      FAPin <- fapPH1ARMA(cc = cc, n = n, order = order, phiVec = phiVec, thetaVec = thetaVec,
			  case = case, method = method, nsim = nsim2, burnIn = burnIn, simType = simType)

    }

    cat('FAPin:', FAPin, 'and cc:', cc, '\n')
    FAP0 - FAPin

  }

  if (is.matrix(phiVec) | is.matrix(thetaVec)) {
    nsim1 <- max(dim(phiVec)[1], dim(thetaVec)[1])
  } else {
    nsim1 <- 0
  }

  ##cat('phiVec:', phiVec, 'thetaVec:', thetaVec, sep = ', ')

  uniroot(root.finding, interval, FAP0 = FAP0, n = n, order = order, phiVec = phiVec, thetaVec = thetaVec, case = case, method = method,
          nsim1 = nsim1, nsim2 = nsim, burnIn = burnIn, simType = simType, seed = seed)$root




}


getCCPH1ARMASimMissingValue <- function(FAP0 = 0.1, interval = c(1, 4), n = 50, mvVec = NULL, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL,
                                        case = 'U', method = 'Method 3', nsim = 1000, burnIn = 50, simType = 'Matrix', logliktol = 1e-2, seed = 12345) {

  root.finding <- function(FAP0, cc, n, mvVec, order, phiVec, thetaVec, case, method, nsim1, nsim2, burnIn, simType, logliktol, seed) {

    set.seed(seed)

    if (nsim1 > 0) {

      FAPin <- lapply(1:nsim1, function(X) {

        phiVecTmp <- phiVec[X, ]
        if (length(phiVecTmp) == 0) phiVecTmp <- rep(0, order[1])

        thetaVecTmp <- thetaVec[X, ]
        if (length(thetaVecTmp) == 0) thetaVecTmp <- rep(0, order[3])

        #cat('AR:', phiVecTmp, '\n')
        #cat('checkAR:', all(abs(outAR[i, ]) < 1), '\n')

        #cat('MA:', thetaVecTmp, '\n')
        #cat('checkMA:', all(abs(outMA[i, ]) < 1), '\n')

        fapPH1ARMAMissingValue(cc = cc, n = n, mvVec = mvVec, order = order, phiVec = phiVecTmp, thetaVec = thetaVecTmp,
                   case = case, method = method, nsim = nsim2, burnIn = burnIn, simType = simType, tol = logliktol)

      })

      FAPin <- mean(unlist(FAPin))

    } else {


      FAPin <- fapPH1ARMAMissingValue(cc = cc, n = n, mvVec = mvVec, order = order, phiVec = phiVec, thetaVec = thetaVec,
                          case = case, method = method, nsim = nsim2, burnIn = burnIn, simType = simType, tol = logliktol)

    }

    cat('FAPin:', FAPin, 'and cc:', cc, '\n')
    FAP0 - FAPin

  }

  if (is.matrix(phiVec) | is.matrix(thetaVec)) {
    nsim1 <- max(dim(phiVec)[1], dim(thetaVec)[1])
  } else {
    nsim1 <- 0
  }

  ##cat('phiVec:', phiVec, 'thetaVec:', thetaVec, sep = ', ')

  uniroot(root.finding, interval, FAP0 = FAP0, n = n, mvVec = mvVec, order = order, phiVec = phiVec, thetaVec = thetaVec, case = case, method = method,
          nsim1 = nsim1, nsim2 = nsim, burnIn = burnIn, simType = simType, logliktol = logliktol, seed = seed)$root




}

getCC <- function(FAP0 = 0.1, interval = c(1, 4), n = 50, mvVec = NULL, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, case = 'U',
              method = 'Method 3', nsimCoefs = 100, nsimProcess = 1000, burnIn = 50, simType = 'Matrix', logliktol = 1e-2, seed = 12345) {


  if (case == 'K') {
    if (simType == 'Matrix') {
      set.seed(seed)
      sigMat <- SigmaMat(n, order = order, phiVec = phiVec, thetaVec = thetaVec, sigma2 = 1, burnIn = burnIn)
      out <- qmvnorm(1 - FAP0, tail = 'both.tails', sigma = sigMat$SigmaMat / sigMat$gamma0)$quantile
    } else {
      if (is.null(mvVec)) {
        out <- getCCPH1ARMASim(FAP0, interval, n, order, phiVec = phiVec, thetaVec = thetaVec, case = case, method = method,
                               nsim = nsimProcess, burnIn = burnIn, simType = simType, seed = seed)
      } else {
        out <- getCCPH1ARMASimMissingValue(FAP0, interval, n, mvVec, order, phiVec = phiVec, thetaVec = thetaVec, case = case, method = method,
                               nsim = nsimProcess, burnIn = burnIn, simType = simType, logliktol = logliktol, seed = seed)
      }

    }
  } else if (case == 'U') {

    if (is.null(mvVec)) {
      CoefDist <- simCoefDist(n, order, phiVec, thetaVec, method, nsim = nsimCoefs, burnIn = burnIn, seed = seed, simType = simType)

      out <- getCCPH1ARMASim(FAP0, interval, n, order, phiVec = CoefDist$phiVec, thetaVec = CoefDist$thetaVec, case = case, method = method,
                             nsim = nsimProcess, burnIn = burnIn, simType = simType, seed = seed)
    } else {

      CoefDist <- simCoefDistMissingValue(n, mvVec, order, phiVec, thetaVec, method, nsim = nsimCoefs, burnIn = burnIn, tol = logliktol, seed = seed, simType = simType)

      out <- getCCPH1ARMASimMissingValue(FAP0, interval, n, mvVec, order, phiVec = CoefDist$phiVec, thetaVec = CoefDist$thetaVec, case = case, method = method,
                             nsim = nsimProcess, burnIn = burnIn, simType = simType, logliktol = logliktol, seed = seed)
    }

  }

  out

}


PH1ARMA <- function(X, cc = NULL, FAP0 = 0.1, order = NULL, plot.option = TRUE, interval = c(1, 4),
                    case = 'U', method = 'Method 3', nsimCoefs = 100, nsimProcess = 1000, burnIn = 50, simType = 'Matrix',
                    iteractive = TRUE, ccReplace = TRUE, logliktol = 1e-2, seed = 12345) {

  Y <- X

  n <- length(Y)

  if (is.null(cc)) {

    if (is.null(order)) {

      if (method == 'Method 1' | method == 'Method 3') {
        model <- auto.arima(Y, method = 'CSS-ML')
      } else if (method == 'Method 2') {
        model <- auto.arima(Y, method = 'CSS')
      }

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

    cc <- getCC(FAP0 = FAP0, interval = interval, n, order = order, phiVec = phiVec, thetaVec = thetaVec, case = case,
      method = method, nsimCoefs = nsimCoefs, nsimProcess = nsimProcess, burnIn = burnIn, simType = simType, seed = seed)

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

