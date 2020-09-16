# Overview
The purpose of PH1ARMA is to build Phase I individual control chart with an ARMA model.  

# Installation

## Install from GitHub

PH1ARMA is still under development, so we recommend users installing it through Github as follows

``` r
install.packages("devtools")
devtools::install_github("bolus123/PH1ARMA")
```

Note that for Windows users,  Rtools may need to be installed in advance.  Please choose the right version of Rtools which is corresponding to your R and Rstudio.  The detailed instruction is introduced: https://cran.r-project.org/bin/windows/Rtools/

For Mac and Linux users, please follow the instruction: https://www.r-project.org/nosvn/pandoc/devtools.html

## Install from local

Users can also download our release, PH1ARMA_x.y.z.tar.gz, from our homepage on Github and then install it from your local path as follows
``` r
install.packages('path_to_file/PH1ARMA_x.y.z.tar.gz', repos = NULL, type="source")
```


# Usage

Before using any functions, PH1ARMA may need to be loaded into R

``` r
library(PH1ARMA)
```

PH1ARMA provides a function to build Phase I individual chart with an ARMA model as follows

``` r
data(chambers_data)
X <- chambers_data ^ (1/3)
PH1ARMA(X)
```

Also, PH1ARMA provides a function to get the corrected charting constant as follows

``` r
# double simulation gets involved if estimators are unknown
getCC(FAP0 = 0.1, double.sim = TRUE)

# single simulation gets involved if estimators are known
getCC(FAP0 = 0.1, double.sim = FALSE)
```

More details are on the manual.
