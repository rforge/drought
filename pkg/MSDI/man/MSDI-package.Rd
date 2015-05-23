\name{MSDI-package}
\alias{MSDI-package}
\alias{MSDI}
\docType{package}
\title{
Multivariate Standardized Drought Index (MSDI)
}
\description{
Compute the Multivariate Standardized Drought Index (MSDI)
}
\details{
\tabular{ll}{
Package: \tab MSDI\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-05-23\cr
License: \tab GPL-3\cr
}
~~ An overview of how to use the package, including the most ~~
~~ important functions ~~
}
\author{
Zengchao Hao

Zengchao Hao <z.hao4univ@gmail.com>
}
\references{
Kao, S. C. and R. S. Govindaraju (2010). A copula-based joint deficit index for droughts. Journal of Hydrology,  380(1-2): 121-134.  

Hao, Z.et al. (2014). Global integrated drought monitoring and prediction system. Scientific Data,  1 (and references therein)
}

\examples{
#' X=runif(120, min = 0, max = 100)
#' Y=runif(120, min = 0, max = 100)
#' fit<-MSDI(X,Y,ts=6) # Compute the 6 month MSDI
#' fit$ProbEmp2 #Get the empirical MSDI
}