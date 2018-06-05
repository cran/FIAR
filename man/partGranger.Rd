\name{partGranger}
\alias{partGranger}
\title{Partial Granger causality
}
\description{Compute partial Granger causality of multivariate timeseries.
}
\usage{
partGranger(data, nx = 1, ny = 1, order=1, perm = FALSE, prob=TRUE, bs = 100)
}
\arguments{
  \item{data}{object containing all observations (rows) and variables (columns)
that are being considered. The variables should be ordered as follows: First
the variables that are supposed to granger cause a set of other variables
(>=1). Then the set of variables (>=1) that are Granger caused by the first set
of variables. Finally, a set of variables to condition on(>=1).
}
  \item{nx}{The number of variables (>=1) that Granger cause a set of other
variables (default = 1), conditioned on a third set of variables (>=1).
}
  \item{ny}{The number of variables (>=1) that are Granger caused by the first
nx variables (default = 1), conditioned on a third set of variables (>=1).
}
  \item{order}{Autoregressive order (>=1) of timeseries. Can be computed using
ARorder().
}
  \item{perm}{Logical. If perm = FALSE (default), only the Granger causality
measure is produced. If perm = TRUE, the Granger test is computed and a
permutation test is performed to do inference.
}
\item{prob}{Logical. If TRUE, the F statistic is returned together with the p-value.
}
  \item{bs}{Number of permutation samples. Only works when perm = TRUE. Default=100
}
}
\value{Partial Granger causality measure F1 plus p -value. 
}
\author{Bjorn Roelstraete}
\references{Guo, S., Seth, A.K., Kendrick, K.M., Zhou, C., Feng, J.(2008).
Partial Granger Causality-Eliminating Exogenous Inputs and Latent Variables.
Journal of Neuroscience Methods. 79-93.
}
\examples{
# Example data with 5 regions x, y, z, q, w
head(grangerdata)

# Calculate AR() order of the data
ARorder(grangerdata, max=10)
 
# Compute partial conditional granger causality of region x to regions y 
# and z, conditional on regions q and w
F <- partGranger(grangerdata, nx=1, ny=2, order=3)

# Compute F and permutation H0 distribution
F <- partGranger(grangerdata, nx=1, ny=2, order=3, perm=TRUE,bs=10)
}
\keyword{ ~Granger }
\keyword{ ~parial }
