\name{getCC}
\alias{getCC}
\title{get a corrected charting constant for Phase I individual control chart with an ARMA model}
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
  \item{order}{ order for ARMA model }
  \item{phiVec}{ autoregressive coeficient vector }
  \item{thetaVec}{ moving average coeficient vector }
  \item{method}{ estimation method for the ARMA model. Default is "CSS-ML". Other options are 'ML' and 'CSS'. }
  \item{double.sim}{ double simulation option.  Default is TRUE. }
  \item{nsimCoefs}{ number of simulation for coeficients.  It is functional when double.sim = TRUE. }
  \item{nsimProcess}{ number of simulation for ARMA processes }
  \item{burnIn}{ number of burn-ins.  When burnIn = 0, the ECM gets involved.  When burnIn is large enough, the ACM gets involved. }
  \item{seed}{ random seed }

}
\description{
Obtain a corrected charting constant.  
}
\examples{
%\dontrun{

# get charting constant with AR(1) model with coeficient = 0.5 using single simulation
getCC(FAP0 = 0.1, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, 
	double.sim = FALSE, nsimProcess = 100)

# get charting constant with MA(1) model with coeficient = 0.5 using single simulation
getCC(FAP0 = 0.1, order = c(0, 0, 1), phiVec = NULL, thetaVec = 0.5, 
	double.sim = FALSE, nsimProcess = 100)

# get charting constant with AR(1) model with coeficient = 0.5 using double simulation
getCC(FAP0 = 0.1, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL, 
	double.sim = TRUE, nsimCoefs = 100, nsimProcess = 100)

# get charting constant with MA(1) model with coeficient = 0.5 using double simulation
getCC(FAP0 = 0.1, order = c(0, 0, 1), phiVec = NULL, thetaVec = 0.5, 
	double.sim = TRUE, nsimCoefs = 100, nsimProcess = 100)

# get charting constant with ARMA(1, 1) model with coeficient = 0.5 using double simulation
getCC(FAP0 = 0.1, order = c(1, 0, 1), phiVec = 0.5, thetaVec = 0.5, 
	double.sim = TRUE, nsimCoefs = 100, nsimProcess = 100)

%}
}
