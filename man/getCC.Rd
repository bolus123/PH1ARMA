\name{getCC}
\alias{getCC}
\title{get a corrected charting constant for Phase I individual control chart with an ARMA model}
 \usage{
getCC(
	FAP0 = 0.1
	,interval = c(1, 4)
	,n = 50
	,order = c(1, 0, 0)
	,phiVec = 0.5
	,thetaVec = NULL
	,case = 'U'
	,method = 'Method 3'
	,nsimCoefs = 100
	,nsimProcess = 1000
	,burnIn = 50
	,simType = 'Matrix'
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
  \item{case}{ known or unknown case.  When case = 'U', the parameters are estimated.  When case = 'K', the parameters are not estimated. }
  \item{method}{ estimation method for the ARMA model. Default is "CSS-ML". Other options are 'ML' and 'CSS'. }
  \item{nsimCoefs}{ number of simulation for coeficients.  It is functional when double.sim = TRUE. }
  \item{nsimProcess}{ number of simulation for ARMA processes }
  \item{burnIn}{ number of burn-ins.  When burnIn = 0, the ECM gets involved.  When burnIn is large enough, the ACM gets involved. }
  \item{simType}{ type of simulation.  When simType = 'Matrix', the simulation is generated using matrix computation.  When simType = 'recursive', the simulation is based on a recursion. }
  \item{seed}{ random seed }

}
\description{
Obtain a corrected charting constant.
}
\examples{
%\dontrun{

# get charting constant with AR(1) model with coeficient = 0.5 using known parameters
getCC(FAP0 = 0.1, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL,
	case = 'K', nsimProcess = 100)

# get charting constant with MA(1) model with coeficient = 0.5 using known parameters
getCC(FAP0 = 0.1, order = c(0, 0, 1), phiVec = NULL, thetaVec = 0.5,
	case = 'K', nsimProcess = 100)

# get charting constant with AR(1) model with coeficient = 0.5 using estimates
getCC(FAP0 = 0.1, order = c(1, 0, 0), phiVec = 0.5, thetaVec = NULL,
	case = 'U', nsimCoefs = 100, nsimProcess = 100)

# get charting constant with MA(1) model with coeficient = 0.5 using estimates
getCC(FAP0 = 0.1, order = c(0, 0, 1), phiVec = NULL, thetaVec = 0.5,
	case = 'U', nsimCoefs = 100, nsimProcess = 100)

# get charting constant with ARMA(1, 1) model with coeficient = 0.5 using estimates
getCC(FAP0 = 0.1, order = c(1, 0, 1), phiVec = 0.5, thetaVec = 0.5,
	case = 'U', nsimCoefs = 100, nsimProcess = 100)

%}
}
