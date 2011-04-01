\name{diffGranger}
\alias{diffGranger}
\title{Directed conditional granger causality 
}
\description{Compute the difference conditional granger causality of
multivariate timeseries. 
}
\usage{
diffGranger(data, nx = 1, ny = 1, order=1, boot = FALSE, bs = 100, p= .05)
}
\arguments{
  \item{data}{object containing all observations (rows) and variables (columns)
that are being considered. The variables should be ordered as follows: First
the variables that are supposed to granger cause a set of other variables
(>=1). Then the set of variables (>=1) that are granger caused by the first set
of variables. Finally, a set of variables to condition on(>=1).
}
  \item{nx}{The number of variables (>=1) that granger cause a set of other
variables (default = 1), conditioned on a third set of variables (>=1).
}
  \item{ny}{The number of variables (>=1) that are granger caused by the first
nx variables (default = 1), conditioned on a third set of variables (>=1).
}
  \item{order}{Autoregressive order (>=1) of timeseries. Can be computed using
ARorder().
}
  \item{boot}{Logical. If boot = FALSE (default), only the granger causality
measure is produced. If boot = TRUE, the granger test is performed and a
bootstrap is performed to generate the H0 distribution (based on function
tsboot() from the package boot).
}
  \item{bs}{Number of bootstrap samples (Default=100). Only works when boot =
TRUE.
}
\item{p}{Uncorrected threshold for significance. This threshold gets Bonferroni corrected.
}
}
\details{The total linear dependence between X and Y can be divided in three
components: a directed influence from X to Y, a directed influence from Y to X
and an undirected instantaneous influence between them. The difference granger
causality from X to Y computed in the function diff.granger() subtracts the
conditional granger causality from Y to X from the conditional granger
causality from X to Y. This can be used as a measure of how much stronger
(weaker) one directed influence is compared to the opposite directed influence. 
}
\references{Roebroeck, A., Formisano, E., Goebel, R. (2005). Mapping directed influence over the brain using Granger causality. NeuroImage 230-242.
}
\value{Difference granger causality measure. When boot=TRUE, the bias and standard error of the bootstrap statistic are also printed.}

\author{Bjorn Roelstraete}
\seealso{condGranger, pdiffGranger}
\examples{
# Example data with 5 regions x, y, z, q, w
head(grangerdata)

# Calculate AR() order of the data
ARorder(grangerdata, max=10)
 
# Compute difference conditional granger causality of region x to regions y o
# and z, conditional on regions q and w
F <- diffGranger(grangerdata, nx=1, ny=2, order=3)

# Compute F and bootstrap H0 distribution
F <- diffGranger(grangerdata, nx=1, ny=2, order=3, boot=TRUE)
}
\keyword{ ~Granger }
\keyword{ ~conditional }
