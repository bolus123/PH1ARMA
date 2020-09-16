\name{PH1ARMA}
\alias{PH1ARMA}
\title{Build Phase I individual control chart with an ARMA model using a corrected charting constant}
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
  \item{X}{ input and it must be a vector }
  \item{cc}{ nominal Phase I charting constant. If this is given, the function will not re-compute the charting constant. }
  \item{FAP0}{ nominal false Alarm Probabilty in Phase I }
  \item{order}{ order for ARMA model }
  \item{method}{ estimation method for the ARMA model. Default is "CSS-ML". Other options are 'ML' and 'CSS'. }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{interval}{ searching range of charting constants for the exact method }
  \item{double.sim}{ double simulation option.  Default is TRUE. }
  \item{nsimCoefs}{ number of simulation for coeficients.  It is functional when double.sim = TRUE. }
  \item{nsimProcess}{ number of simulation for ARMA processes }
  \item{burnIn}{ number of burn-ins.  When burnIn = 0, the ECM gets involved.  When burnIn is large enough, the ACM gets involved. }
  \item{seed}{ random seed }
}
\description{
Build a Phase I individual control chart with an ARMA model. The charting constant is correted by this approach.
}
\examples{
%\dontrun{

# load the data in the package as an example
data(data)
X <- data ^ (1/3)

# An example using the default setting whose FAP0 = 0.1
PH1ARMA(X, nsimCoefs = 100, nsimProcess = 100)

# An example using the default setting for FAP0 = 0.05
PH1ARMA(X, FAP0 = 0.05, nsimCoefs = 100, nsimProcess = 100)

# An example using the singular simualtion with FAP0 = 0.1
PH1ARMA(X, double.sim = FALSE, nsimProcess = 100)

# An example using the singular simualtion with FAP0 = 0.05
PH1ARMA(X, FAP0 = 0.05, double.sim = FALSE, nsimProcess = 100)

%}
}
