\name{getCC}
\alias{getCC}
\title{get Phase I corrected charting constant}
 \usage{
getCC(
	FAP0 = 0.1
	,interval = c(1, 4)
	,n = 30
	,order = c(1, 0, 0)
	,phiVec = 0.5
	,thetaVec = NULL
	,method = 'CSS-ML'
	,double.sim = TRUE
	,nsimCoefs = 1000
	,nsimProcess = 1000
	,burnIn = 1000
	,seed = 12345
)
}
\arguments{
  \item{FAP0}{ nominal false Alarm Probabilty in Phase 1 }
  \item{interval}{ searching range of charting constants for the exact method }
  \item{n}{ number of observations }
  \item{order}{ order of ARMA model }
  \item{phiVec}{ phi vector }
  \item{thetaVec}{ theta vector }
  \item{method}{ estimation method }
  \item{double.sim}{ double simulation }
  \item{nsimCoefs}{ number of simulation for coeficients }
  \item{nsimProcess}{ number of simulation for processes with ARMA }
  \item{burnIn}{ number of burn-ins }
  \item{seed}{ random seed }

}
\description{
Obtain a corrected charting constant.  
}
\examples{
%\dontrun{

# get charting constant with AR(1) model with coeficient = 0.5 using single simulation
getCC(FAP0 = 0.1, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, double.sim = FALSE, nsimProcess = 100)

# get charting constant with MA(1) model with coeficient = 0.5 using single simulation
getCC(FAP0 = 0.1, order = c(0, 0, 1), phiVec = NULL, thetaVec = 0.5, double.sim = FALSE, nsimProcess = 100)

# get charting constant with AR(1) model with coeficient = 0.5 using double simulation
getCC(FAP0 = 0.1, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, double.sim = TRUE, nsimCoefs = 100, nsimProcess = 100)

# get charting constant with MA(1) model with coeficient = 0.5 using double simulation
getCC(FAP0 = 0.1, order = c(0, 0, 1), phiVec = NULL, thetaVec = 0.5, double.sim = TRUE, nsimCoefs = 100, nsimProcess = 100)

# get charting constant with ARMA(1, 1) model with coeficient = 0.5 using double simulation
getCC(FAP0 = 0.1, order = c(1, 0, 1), phiVec = 0.5, thetaVec = 0.5, double.sim = TRUE, nsimCoefs = 100, nsimProcess = 100)

%}
}
