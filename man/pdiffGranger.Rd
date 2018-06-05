\name{pdiffGranger}
\alias{pdiffGranger}
\title{partial difference Granger causality
}
\description{Compute partial difference conditional Granger causality of
multivariate timeseries.
}
\usage{
pdiffGranger(data, nx = 1, ny = 1, order=1, perm = FALSE, bs = 100)
}
\arguments{
  \item{data}{object containing all observations (rows) and variables (columns)
that are being considered. The variables should be ordered as follows: First
the variables that are supposed to Granger cause a set of other variables
(>=1). Then the set of variables (>=1) that are Granger caused by the first set
of variables. Finally, a set of variables to condition on(>=1).
}
  \item{nx}{The number of variables (>=1) that are supposed to Granger cause a
set of other variables (default = 1), conditioned on a third set of variables
(>=1).
}
  \item{ny}{The number of variables (>=1) that are supposed to be Granger
caused by the first nx variables (default = 1), conditioned on a third set of
variables (>=1).
}
  \item{order}{Autoregressive order (>=1) of timeseries. Can be computed using
ARorder().
}
  \item{perm}{Logical. If perm = FALSE (default), only the Granger causality
measure is produced. If perm = TRUE, the Granger test is computed and a
permutation is performed to generate the H0 distribution.
}
  \item{bs}{Number of permutation samples. Default=100
}
}
\details{The total linear dependence between X and Y can be divided in three
components: a directed influence from X to Y, a directed influence from Y to X
and an undirected instantaneous influence between them. The difference Granger
causality from X to Y computed in the function diff.Granger() subtracts the
partial conditional Granger causality from Y to X from the partial conditional
Granger causality from X to Y. This can be used as a measure of how much
stronger (weaker) one directed influence is compared to the opposite directed
influence. 
}
\value{Partial difference Granger causality measure and p value.}
\seealso{diffGranger}
\author{Bjorn Roelstraete}
\examples{
# Example data with 5 regions x, y, z, q, w
head(grangerdata)

# Calculate AR() order of the data
ARorder(grangerdata, max=10)
 
# Compute partial difference conditional Granger causality of region x to 
# regions y and z, conditional on regions q and w
F <- pdiffGranger(grangerdata, nx=1, ny=2, order=3)

# Compute F and permutation H0 distribution
F <- pdiffGranger(grangerdata, nx=1, ny=2, order=3, perm=TRUE, bs=50)
}
\keyword{ ~Granger }
\keyword{ ~partial }
