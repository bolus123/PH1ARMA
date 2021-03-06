\name{PH1ARMA}
\alias{PH1ARMA}
\title{Build Phase I individual control chart with an ARMA model using a corrected charting constant}
\usage{
PH1ARMA(
	X
	,cc = NULL
	,FAP0 = 0.1
	,order = NULL
	,method = 'Method 3'
	,plot.option = TRUE
	,interval = c(1, 4)
	,case = 'U'
	,nsimCoefs = 100
	,nsimProcess = 1000
	,burnIn = 50
	,simType = 'Matrix'
	,seed = 12345
)
}
\arguments{
  \item{X}{ input and it must be a vector }
  \item{cc}{ nominal Phase I charting constant. If this is given, the function will not re-compute the charting constant. }
  \item{FAP0}{ nominal false Alarm Probabilty in Phase I }
  \item{order}{ order for ARMA model }
  \item{method}{ estimation method for the control chart. When method = 'Method 3' is maximum likehood estimations plus method of moments. Other options are 'Method 1' which is pure MLE and 'Method 2' which is pure CSS. }
  \item{plot.option}{TRUE - draw a plot for the process; FALSE - Not draw a plot for the process}
  \item{interval}{ searching range of charting constants for the exact method }
  \item{case}{ known or unknown case.  When case = 'U', the parameters are estimated }
  \item{nsimCoefs}{ number of simulation for coeficients.  It is functional when double.sim = TRUE. }
  \item{nsimProcess}{ number of simulation for ARMA processes }
  \item{burnIn}{ number of burn-ins.  When burnIn = 0, the ECM gets involved.  When burnIn is large enough, the ACM gets involved. }
  \item{simType}{ type of simulation.  When simType = 'Matrix', the simulation is generated using matrix computation.  When simType = 'Recursive', the simulation is based on a recursion. }
  \item{seed}{ random seed }
}
\description{
Build a Phase I individual control chart with an ARMA model. The charting constant is correted by this approach.
}
\examples{
%\dontrun{

# load the data in the package as an example
data(chambers_data)
X <- chambers_data ^ (1/3)

# An example using the default setting whose FAP0 = 0.1
PH1ARMA(X, nsimCoefs = 100, nsimProcess = 100)

# An example using the default setting for FAP0 = 0.05
PH1ARMA(X, FAP0 = 0.05, nsimCoefs = 100, nsimProcess = 100)

# An example using known parameters with FAP0 = 0.1
PH1ARMA(X, case = 'K', nsimProcess = 100)

# An example using the single simualtion with FAP0 = 0.05
PH1ARMA(X, FAP0 = 0.05, case = 'K', nsimProcess = 100)

# When users get an error message about the size of matrix,
# the function needs to use the alternative simulation type as follows
PH1ARMA(X, FAP0 = 0.05, nsimCoefs = 100, nsimProcess = 100, simType = 'Recursive')


%}
}
