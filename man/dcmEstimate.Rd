\name{dcmEstimate}
\alias{dcmEstimate}
\title{DCM estimator
}
\description{Estimate parameters of a bilinear DCM}
\usage{
dcmEstimate(DCM,ts)
}
\arguments{
 \item{DCM}{DCM object.}
 \item{ts}{Timeseries to fit the model to.}
 }
\value{Returns posterior parameter values:
  \item{DCM$A}{posterior anatomical connections from column j to row i.}
  \item{DCM$B}{posterior functional connections from column j to row i for
every input [,,k].}
  \item{DCM$C}{posterior input connections from input k (row) to region l
(column).}
  \item{DCM$H}{posterior hemodynamic parameters.}
  \item{DCM$Cp}{posterior parameter covariance.}
  \item{DCM$Ce}{posterior error covariance.}
}

\references{Friston et al. 2003. Dynamic Causal Modeling. Neuroimage, 19,
1273-1302.
}
\author{Bjorn Roelstraete}
\seealso{dcmGenerate}
\examples{
# Estimate posterior parameters of model DCMex with data DCMex$sim
 ts <-  dcmGenerate(DCMex, SNR=1, ar=.2, names=c('V1','V2','V3'))
 DCMex <- dcmEstimate(DCMex, ts)

#Posterior anatomical connections
 DCMex$A

#Posterior functional connections
 DCMex$B

#Posterior input connections
 DCMex$C

}
\keyword{ ~DCM }
