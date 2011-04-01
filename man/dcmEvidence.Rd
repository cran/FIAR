\name{dcmEvidence}
\alias{dcmEvidence}
\title{modelfit}
\description{Compute AIC and BIC of a DCM
}
\usage{
dcmEvidence(DCM, ts)
}
\arguments{
 \item{DCM}{DCM object.}
 \item{ts}{Timeseries to fit the model to.}
 }

\value{Creates 2 extra fields DCM$AIC and DCM$BIC}
\seealso{spm.dcm.estimate}
\author{Bjorn Roelstraete}
\examples{
# Compute how well the model DCMex fits the timeseries DCMex$sim
 ts <-  dcmGenerate(DCMex, SNR=1, ar=.2, names=c('V1','V2','V3'))
 DCMex <- dcmEvidence(DCMex, ts)

 DCMex$AIC
 DCMex$BIC


}
\keyword{ ~DCM }
