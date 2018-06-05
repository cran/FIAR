\name{dcmGenerate}
\alias{dcmGenerate}
\title{DCM timeseries generator
}
\description{Generate simulated timeseries from a specified DCM
}
\usage{dcmGenerate(DCM, SNR = 0, ar = 0,names=DCM$names)
}
\arguments{
  \item{DCM}{A DCM list containing all model and experimental parameters. This
list can be constructed using dcm_param or manually.
}
  \item{SNR}{SNR of the timeseries. The number represents sd(signal)/sd(noise).
If SNR=0 the pure signal is generated. 
}
  \item{ar}{Autoregressioncoefficient of the noise added. 0 (default) means
white, gaussian noise.
}
\item{names}{The names of the variables.}
}
\value{Function creates a field DCM$sim, which contains the simulated
timeseries per timepoint (rows) and region (columns).
}
\author{Bjorn Roelstraete}
\seealso{dcm_param}
\examples{
# Use example DCMex to generate three timeseries V1, V2, V3 with a SNR of 1 and AR(.2)

 ts <-  dcmGenerate(DCMex, SNR=1, ar=.2, names=c('V1','V2','V3'))
 plot(ts[,1], t='b')

}
\keyword{ ~DCM }
