% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CH06_LM.R
\name{FCI_CH06}
\alias{FCI_CH06}
\title{Residual-based test for fractional cointegration (Chen, Hurvich (2006))}
\usage{
FCI_CH06(X, m_peri, m, alpha = 0.05, diff_param = 1)
}
\arguments{
\item{X}{data matrix.}

\item{m_peri}{fixed positive integer for averaging the periodogram, where \code{m_peri>(nbr of series + 3)}}

\item{m}{bandwith parameter specifying the number of Fourier frequencies
used for the estimation, usually \code{floor(1+T^delta)}, where 0<delta<1.}

\item{alpha}{desired significance level. Default is \code{alpha=0.05}.}

\item{diff_param}{integer specifying the order of differentiation in order to ensure stationarity of data, where diff_param-1 are the number of differences. 
Default is \code{diff_param=1} for no differences.}
}
\description{
\code{FCI_CH06} Semiparametric residual-based test for fractional cointegration by Chen, Hurvich (2003).
 Returns test statistic, critical value and testing decision. Null hypothesis: no fractional cointegration.
}
\examples{
T<-1000
series<-FI.sim(T=T, q=2, rho=0.4, d=c(0.1,0.4), B=rbind(c(1,-1),c(0,1)))
FCI_CH06(series, diff_param=1, m_peri=25, m=floor(T^0.65))
series<-FI.sim(T=T, q=2, rho=0.4, d=c(0.4,0.4))
FCI_CH06(series, diff_param=1, m_peri=25, m=floor(T^0.65))
}
\references{
Chen, W. W. and Hurvich, C. M. (2006): Semiparametric estimation of fractional
cointegrating subspaces. The Annals of Statistics, Vol. 34, No. 6, pp. 2939 - 2979.
}
\author{
Christian Leschinski
}
