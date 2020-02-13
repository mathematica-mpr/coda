\name{effectiveSize}
\alias{effectiveSize}
\title{Effective sample size for estimating the mean}
\description{
Sample size adjusted for autocorrelation.
}
\usage{
effectiveSize(x, crosschain=FALSE) 
}
\arguments{
\item{x}{An mcmc or mcmc.list object.}
\item{crosschain}{A logical flag indicating whether to incorporate information from cross-chain variance. Only available when \code{x} is an mcmc.list object with multiple chains.}
}
\details{
For a time series \code{x} of length \code{N}, the standard error of the
mean is the square root of \code{var(x)/n} where \code{n} is the
effective sample size.  \code{n = N} only when there is no
autocorrelation.

Estimation of the effective sample size requires estimating the
spectral density at frequency zero.  This is done by the function
\code{spectrum0.ar}

For a \code{mcmc.list} object, the effective sizes are summed across
chains.  To get the size for each chain individually use
\code{lapply(x,effectiveSize)}.  

When \code{crosschain} is set to \code{TRUE}, spectral densities are not calculated, and instead 
autocorrelations are calculated using a fast Fourier transform, and are adjusted for the estimate of
cross-chain variance.
}
\value{
A vector giving the effective sample size for each column of \code{x}.
}
\seealso{
   \code{\link{spectrum0.ar}}.
}
\references{
For calculations used by the \code{crosschain} option: The Stan Development Team \emph{Stan Reference Manual}, Section 16.4 \href{https://mc-stan.org/docs/2_21/reference-manual/effective-sample-size-section.html}{here}.
}
\keyword{ts}
