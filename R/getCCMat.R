#library(mvtnorm)
#library(pracma)
##############################################
"InvertQ" <- 
function(coef){
 stopifnot(class(coef)=="numeric"||class(coef)=="matrix"||(class(coef)=="array" && (dim(coef)[1]==dim(coef)[2])))
     if (class(coef) == "numeric")
       coef <- array(coef,dim=c(1,1,length(coef)))
     if (class(coef) == "matrix")
       coef <- array(coef,dim=c(NROW(coef),NROW(coef),1))
     k <- dim(coef)[1]
     order <- dim(coef)[3]
     if (order==1)
       ans <- eigen(coef[,,1], symmetric=FALSE, only.values =TRUE)$value 
     else{
        blockMat <- matrix(numeric((k*order)^2),k*order,k*order)
        blockMat[1:k,] <- coef
        Imat <- diag((order-1)*k)
        blockMat[((k+1):(k*order)),1:((order-1)*k)] <- Imat
        ans <- eigen(blockMat, symmetric=FALSE, only.values =TRUE)$value 
     }
   MaxEigenvalue <- max(Mod(ans))
   if (MaxEigenvalue >= 1) 
     return( warning("check stationary/invertibility condition !"))
}
##############################################

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

parsMat <- function(n, parsVec, norder = 1) {
    Check <- InvertQ(parsVec)
    if (class(Check) == 'character') {
        NULL
    } else if (is.null(Check)) {
        Mat <- diag(n)
        for (i in 1:norder) {
            Mat <- Mat + Diag(rep(parsVec[i], n - i), k = -i)
        }
        Mat
    }
}

SigmaMat <- function(n, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, sigma2 = 1, burnIn = 30) {
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
      out <- out[-c(1:burnIn), -c(1:burnIn)] / out[n + burnIn, n + burnIn]
    }

    list(SigmaMat = out, sqrtSigmaMat = sqrtm(out)$B, gamma0 = gamma0)
}

simARMAProcessMat <- function(n, sqrtSigMat, innovDist = 'norm', innovPars = c(0, 1)) {

	sqrtSigMat %*% matrix(simInnov(n, sigma = 1, XSim = innovDist, XPars = innovPars), ncol = 1)
	
}

getCCPH1ARMAMatKnown <- function(FAP0 = 0.1, n = 30, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL,
                         MatType = 'ECM', seed = 12345) {
    set.seed(seed)
    if (MatType == 'ECM') {
        sigMat <- SigmaMat(n, order = order, phiVec = phiVec, thetaVec = thetaVec, sigma2 = 1, burnIn = 0)$SigmaMat
    } else if (MatType == 'ACM') {
        sigMat <- SigmaMat(n, order = order, phiVec = phiVec, thetaVec = thetaVec, sigma2 = 1, burnIn = 30)$SigmaMat
    }
    qmvnorm(1 - FAP0, tail = 'both.tails', sigma = sigMat)$quantile
}

simCoefDistMat <- function(n, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, method = 'Method 3',
                        MatType = 'ECM', nsim = 1000, seed = 12345) {
  
  set.seed(seed)
  
  outAR <- matrix(NA, nrow = nsim, ncol = length(phiVec))
  outMA <- matrix(NA, nrow = nsim, ncol = length(thetaVec))
  outMean <- rep(NA, nsim)
  outGamma <- rep(NA, nsim)
  
  if (MatType == 'ACM') {
    sigMat <- SigmaMat(n, order, phiVec, thetaVec, sigma2 = 1, burnIn = 30)
  } else if (MatType == 'ECM') {
    sigMat <- SigmaMat(n, order, phiVec, thetaVec, sigma2 = 1, burnIn = 0)
  }

  for (i in 1:nsim) {
    
      flg <- 1
    
      while (flg == 1) {
      
        sim <- simARMAProcessMat(n, sigMat$sqrtSigmaMat, innovDist = 'norm', innovPars = c(0, 1))
        if (method == 'Method 1' | method == 'Method 3') {
            model <- try(arima(sim, order = order, method = 'CSS-ML'), silent = TRUE)
        } else if (method == 'Method 2') {
            model <- try(arima(sim, order = order, method = 'CSS'), silent = TRUE)
        }
        
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


fapPH1ARMAMatUnknown <- function(cc = 3, n = 30, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, 
  MatType = 'ECM', method = 'Method 3', nsim = 1000) {
					   
	if (MatType == 'ACM') {
    sigMat <- SigmaMat(n, order, phiVec, thetaVec, sigma2 = 1, burnIn = 30)
  } else if (MatType == 'ECM') {
    sigMat <- SigmaMat(n, order, phiVec, thetaVec, sigma2 = 1, burnIn = 0)
  }			   
				   
	out <- lapply(1:nsim, function(X){
		flg <- 1
		while(flg == 1) {
			sim <- simARMAProcessMat(n, sigMat$sqrtSigmaMat, innovDist = 'norm', innovPars = c(0, 1))
				if (method == 'Method 1') {
					m1 <- try(arima(sim, order = order, method = 'CSS-ML'), silent = TRUE)
					if (class(m1) == 'try-error') {
						flg = 1
					} else {

            if (!is.null(phiVec)) {
              phiVec1 <- m1$coef[1:order[1]]
            } else {
              phiVec1 <- NULL
            }
            
            if (!is.null(thetaVec)) {
              thetaVec1 <- m1$coef[(order[1] + 1):(order[1] + order[3])]
            } else {
              thetaVec1 <- NULL
            }

            intercept <- m1$coef[length(m1$coef)]
            mu0 <- intercept * (1 - sum(phiVec1))

            sigMat1 <- SigmaMat(n, order, phiVec1, thetaVec1, sigma2 = m1$sigma2, burnIn = 30)

						sim <- (sim - mu0) / sqrt(sigMat1$gamma0)
						flg <- 0
					}
				} else if (method == 'Method 2' | method == 'Method 3') {
					sim <- (sim - mean(sim)) / sd(sim)
					flg <- 0
				}
		}
		sum(-cc <= sim & sim <= cc) != n
	})
	
	mean(unlist(out))

}

getCCPH1ARMAMatUnknown <- function(FAP0 = 0.1, interval = c(1, 4), n = 30, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, 
                        MatType = 'ECM', method = 'Method 3', double.sim = TRUE, nsimCoefs = 1000, nsimProcess = 1000, seed = 12345) {

  root.finding <- function(cc, FAP0 = 0.1, n = 30, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, 
                        MatType = 'ECM', method = 'Method 3', nsimCoefs = 1000, nsimProcess = 1000, seed = 12345) {

      set.seed(seed)

      FAPin <- mean(unlist(lapply(1:nsimCoefs, function(X) {
        fapPH1ARMAMatUnknown(cc, n, order = order, phiVec = phiVec[X, ], thetaVec = thetaVec[X, ], 
           MatType = MatType, method = method, nsim = nsimProcess)}
      )))
     
      cat('FAPin:', FAPin, 'cc:', cc, '\n')
      FAP0 - FAPin

  }

  if (double.sim == TRUE) {
    CoefDist <- simCoefDistMat(n, order = order, phiVec = phiVec, thetaVec = thetaVec, method = method,
                        MatType = MatType, nsim = nsimCoefs, seed = seed)
  } else {
    nsimCoefs = 1
  }

  uniroot(root.finding, interval, FAP0 = FAP0, n = n, order = order, phiVec = CoefDist$phiVec, thetaVec = CoefDist$thetaVec, 
          MatType = MatType, method = method, nsimCoefs = nsimCoefs, nsimProcess = nsimProcess, seed = seed)$root

}

getCCPH1ARMAMatUnknown(FAP0 = 0.1, interval = c(1, 4), n = 50, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, 
                        MatType = 'ECM', method = 'Method 1', double.sim = TRUE, nsimCoefs = 1000, nsimProcess = 100, seed = 12345)
