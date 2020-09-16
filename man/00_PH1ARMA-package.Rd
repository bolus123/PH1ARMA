\name{PH1ARMA-package}
\alias{00_PH1ARMA-package}
\docType{package}
\title{
Phase I Individual Control Chart with an ARMA model
}
\description{
The utility of this package is in building a Phase I individual control chart with an ARMA model. The Phase I chart is calculated using a Monte Carlo approach.
}
\details{
\tabular{ll}{
Package: \tab PH1ARMA\cr
Type: \tab Package\cr
% Version: \tab 1.0\cr
Date: \tab 2020-09-15\cr
License: \tab GPL (>= 3)\cr
}
Function 'PH1ARMA' builds a Phase I individual control chart with an ARMA model using a correct charting constant.
Function 'getCC' gets a charting constant for the data with an ARMA model.
}
\author{ Yuhui Yao, Subha Chakraborti, Xin Thomas Yang, Jason Parton, Dwight Lewis.

Maintainer: Yuhui Yao <yyao17@cba.ua.edu>
}
\references{
abc
}
\examples{
%\dontrun{

# Build a Phase I individual control chart with an ARMA model
data(chambers_data)
X <- chambers_data ^ (1/3)
PH1ARMA(X, nsimCoefs = 100, nsimProcess = 100)

%}
}
