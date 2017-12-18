\name{invClustEM}
\alias{invClustEM}
\alias{EM}
\alias{EM1D}
\alias{getQuality2D}
\alias{getQuality}
\alias{selectSNPsLD}
\alias{invFreq}
\alias{invGenotypes}

\title{
Internal functions of \code{invClust}
}
\description{
Internal functions and Expectation Maximization routines not intended for user use.
}
\usage{
EM(y, mu10, mu20, sigma0, p0, q0, x, muGroup0, sigmaGroup0, piGroup0, tol = 10^(-5), it = 1000)
EM1D(y, mu10, mu20, sigma0, p0, q0, x, muGroup0, sigmaGroup0, piGroup0, tol = 10^(-5), it = 1000)
}
\arguments{
  \item{y}{MDS compoents}
  \item{mu10}{initial condition for mean of first homozygous group}
  \item{mu20}{initial condition for mean of second homozygous group}
  \item{sigma0}{initial condition for covariant matrix between MDS components}
  \item{p0}{initial condition for first haplotype allele}
  \item{q0}{initial condition for second haplotype allele}
  \item{x}{first genome-wide PCA component}
  \item{muGroup0}{initial condition for geographical subpopulation groups}
  \item{sigmaGroup0}{initial condition for standard deviation of groups}
  \item{piGroup0}{initial condition for probability of belonging to each group}
  \item{tol}{tolerace of convergence}
  \item{it}{maxium number of iterations}
}

\seealso{
invClust
}


\keyword{interanal}
\keyword{invClust}