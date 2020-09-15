\name{PH1ARMA-package}
\alias{00_PH1ARMA-package}
\docType{package}
\title{
Phase I Individual Control Chart with an ARMA model
}
\description{
The utility of this package is in building a Shewhart-type control chart based on new methods for the charting statistics, subgroup means. The Phase I chart is based on multivariate normal/t distribution.
}
\details{
\tabular{ll}{
Package: \tab PH1ARMA\cr
Type: \tab Package\cr
% Version: \tab 1.0\cr
Date: \tab 2018-01-05\cr
License: \tab GPL (>= 3)\cr
}
Function 'PH1XBAR' builds a Phase I Shewhart X-bar control chart with a correct charting constant.
Function 'getCC' gets a charting constant for the data with a balanced one-way random effects model.
}
\author{ Yuhui Yao, Subha Chakraborti, Xin Thomas Yang, Jason Parton, Dwight Lewis.

Maintainer: Yuhui Yao <yyao17@cba.ua.edu>
}
\references{
abc
}
\examples{
%\dontrun{

# Build a Phase I basic Shewhart control chart
data(data)
X <- data ^ (1/3)
PH1ARMA(X)


%}
}
