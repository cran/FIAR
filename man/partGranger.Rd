\name{partGranger}
\alias{partGranger}
\title{Partial Granger causality
}
\description{Compute partial Granger causality of multivariate timeseries.
}
\usage{
partGranger(data, nx = 1, ny = 1, order=1, boot = FALSE, bs = 100,b=NA,p=.05)
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
measure is produced. If boot = TRUE, the granger test is computed and a
bootstrap is performed to do inference.
}
  \item{bs}{Number of bootstrap samples. Only works when p = TRUE. Default=100
}
\item{b}{Windowsize. Only works when boot = TRUE.
}
\item{p}{Uncorrected threshold for significance. This threshold gets Bonferroni corrected.
}
}
\value{Partial Granger causality measure F1. When boot=TRUE, the original measurements and a probability matrix is returned. 
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

# Compute F and bootstrap H0 distribution
F <- partGranger(grangerdata, nx=1, ny=2, order=3, boot=TRUE,bs=10,b=10)
}
\keyword{ ~Granger }
\keyword{ ~parial }
