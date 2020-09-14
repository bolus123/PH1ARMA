library(forecast)

simInnov <- function(n, sigma = 1, XSim = 'norm', XPars = c(0, 1)) {

  if (XSim == "norm") {

    me = XPars[1];
    std = sqrt(XPars[2]);

    pars <- c(me, std)

    rgen <- rnorm

  } else if (XSim == "t") {

    me = 0;
    std = sqrt(XPars[1] / (XPars[1] - 2));

    pars <- c(XPars[1])

    rgen <- rt

  } else if (XSim == "gamma") {

    me = XPars[1] * XPars[2];
    std = sqrt(XPars[1] * XPars[2] ^ 2);

    pars <- c(XPars[1], 1 / XPars[2])

    rgen <- rgamma


  } else if (XSim == "beta") {

    me = XPars[1] / (XPars[1] + XPars[2]);
    std = sqrt(XPars[1] * XPars[2] / ((XPars[1] + XPars[2]) ^ 2) / (XPars[1] + XPars[2] + 1));

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

    } else {
      i <- i - 1
    }

  }

  list(phiVec = outAR, thetaVec = outMA)

}


fapPH1ARMA <- function(cc = 3, n = 30, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL,
                         nsim = 1000, burnIn = 1000) {

  out <- lapply(1:nsim, function(X){
    sim <- simARMAProcess(n, order, phiVec, thetaVec, sigma = 1,
                          innovDist = 'norm', innovPars = c(0, 1), burnIn = burnIn)
    sim <- (sim - mean(sim)) / sd(sim)
    sum(-cc <= sim & sim <= cc) != n
  })

  mean(unlist(out))

}

getCCPH1ARMA <- function(FAP0 = 0.1, interval = c(1, 4), n = 30, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL,
                         nsim1 = 1000, nsim2 = 1000, burnIn = 1000, seed = 12345) {

  root.finding <- function(FAP0, cc, n, order, phiVec, thetaVec, nsim1, nsim2, burnIn, seed) {

    set.seed(seed)

    if (nsim1 > 0) {

      FAPin <- lapply(1:nsim1, function(X) {

        phiVecTmp <- phiVec[X, ]
        if (length(phiVecTmp) == 0) phiVecTmp <- NULL

        thetaVecTmp <- thetaVec[X, ]
        if (length(thetaVecTmp) == 0) thetaVecTmp <- NULL

        fapPH1ARMA(cc, n, order, phiVecTmp, thetaVecTmp, nsim2, burnIn)

      })

      FAPin <- mean(unlist(FAPin))

    } else {

      FAPin <- fapPH1ARMA(cc, n, order, phiVec, thetaVec, nsim2, burnIn)

    }

    cat('FAPin:', FAPin, 'and cc:', cc, '\n')
    FAP0 - FAPin

  }

  uniroot(root.finding, interval, FAP0 = FAP0, n = n, order = order, phiVec = phiVec, thetaVec = thetaVec,
          nsim1 = nsim1, nsim2 = nsim2, burnIn = burnIn, seed = seed)$root


}
