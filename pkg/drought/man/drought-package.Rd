\name{drought-package}
\alias{drought-package}
\alias{drought}
\docType{package}
\title{
 Tools for Drought Modeling and Assessment
}
\description{
Tools to be used for drought monitoring, prediction and analysis.
}
\details{
\tabular{ll}{
Package: \tab drought\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-15\cr
License: \tab GPL-3\cr
}

}
\author{
Zengchao Hao

Maintainer: Zengchao Hao Z.hao4univ@gmail.com
}
\references{

Hao, Z., and V. P. Singh(2015), Drought characterization from a multivariate perspective: A review J. Hydrol., 527

Hao, Z., F. Hao, V. P. Singh, and H. Cheng (2015), Drought: A comprehensive R package for drought monitoring, prediction and analysis, Geophysical Research Abstracts, Vol. 17, EGU2015-8070-3, 2015, EGU General Assembly 2015.

Hao, Z. and V. P. Singh (2012), Entropy based method for drought analysis, Journal of Hydrologic Engineering, 18(7):780 - 786.

Kao, S. C. and R. S. Govindaraju (2010). A copula-based joint deficit index for droughts. Journal of Hydrology,  380(1-2): 121-134.  

Hao, Z.et al. (2014). Global integrated drought monitoring and prediction system. Scientific Data,  1 (and references therein)
}

\examples{
#' X=runif(120, min = 0, max = 100)
#' Y0<-ACCU(X,ts=6) # Compute the 6 month  accumulated series
#' SPI<-SDI(X,ts=6) # Get the standardized drought index (or SPI)
#' Y=runif(120, min = 0, max = 100)
#' fit<-MDI(X,Y,ts=6) # Compute the 6 month multivariate drought index
#' fit$ProbEmp2 #Get the empirical multivariate drought index
}


