\name{PH1ARMA}
\alias{PH1ARMA}
\title{Build Phase I X-bar control chart with a corrected charting constant}
\usage{
PH1ARMA(
	X
	,cc = NULL
	,FAP0 = 0.1
	,order = NULL
	,method = 'CSS-ML'
	,plot.option = TRUE
	,interval = c(1, 4)
	,double.sim = TRUE
	,nsimCoefs = 1000
	,nsimProcess = 1000
	,burnIn = 1000
	,seed = 12345
)
}
\arguments{
  \item{X}{ input and it must be a matrix }
  \item{cc}{ nominal Phase I charting constant. If this is given, the function will not recompute the charting constant. }
  \item{FAP0}{ nominal false Alarm Probabilty in Phase 1 }
  \item{order}{ abc }
  \item{method}{ abc }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{interval}{ searching range of charting constants for the exact method }
  \item{double.sim}{ abc }
  \item{nsimCoefs}{ number of simulation for the exact method }
  \item{nsimProcess}{ number of simulation for the exact method }
  \item{burnIn}{ number of simulation for the exact method }
  \item{seed}{ random seed for the exact method }
}
\description{
Build a Phase I Shewhart X-bar control chart with a variance component model. The charting constant is correted by this approach.
}
\examples{
%\dontrun{

# load the data in the package as an example
data(data)

X <- data ^ (1/3)

# An example using the default setting whose FAP0 = 0.1
#PH1ARMA(X, nsimCoefs = 100, nsimProcess = 100)

# An example for FAP0 = 0.05
#PH1ARMA(X, FAP0 = 0.05, nsimCoefs = 100, nsimProcess = 100)

# An example using the default setting whose FAP0 = 0.1
#PH1ARMA(X, double.sim = FALSE, nsimProcess = 100)

# An example for FAP0 = 0.05
#PH1ARMA(X, FAP0 = 0.05, double.sim = FALSE, nsimProcess = 100)

%}
}
