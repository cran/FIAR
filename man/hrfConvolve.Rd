\name{hrfConvolve}
\alias{hrfConvolve}
\title{Convolution function
}
\description{Convolute a timeseries with any double gamma function.
}
\usage{
hrfConvolve(x = NULL, scans = NA, onsets = c(), durations = c(), 
             rt = NA, SNR = 0, mean = FALSE, a1 = 6, a2 = 12, b1 = 0.9, 
             b2 = 0.9, cc = 0.35)
}
\arguments{
  \item{x}{Single timeseries (default = NULL)
}
  \item{scans}{number of scans
}
  \item{onsets}{onsets of experimental condition
}
  \item{durations}{duration of experimental condition
}
  \item{rt}{repetition time
}
  \item{SNR}{signal to noise ratio of data
}
  \item{mean}{logical if mean is TRUE the timeseries is centered around 0.
}
\item{a1}{parameter of the double gamma function}
\item{a2}{parameter of the double gamma function}
\item{b1}{parameter of the double gamma function}      
\item{b2}{parameter of the double gamma function}
\item{cc}{parameter of the double gamma function}
}
\details{The function is an extension of the fmri.stimulus function in the
`fmri' pacakge (see ref.). If x = NULL, the to be convolved stimulusfunction
can be specified with the parameters 'scans', 'onsets', 'durations', 'rt', and
'SNR'. If x is entered, the timeseries x is convolved and the other parameters
need not be specified. The default convolution function is the canonical HRF,
but can be altered by changing the parameters of the double gamma function.
}
\value{returns convolved timeseries. The timeseries is convolved with a mixture
of 2 gamma functions (default = canonical HRF).
}
\references{Polzehl, J. and Tabelow, K. (2007) _fmri: A Package for Analyzing
     fmri Data_, R News, 7:13-17 .
}
\author{Bjorn Roelstraete}
\examples{
# Specify a stimulusfunction without noise and convolve with canonical HRF

hrfConvolve(scans = 240, onsets = c(0,60,120,180), durations = c(30), 
             rt = 3, SNR = 0)

# Convolve a (part of a) timeseries with a canonical HRF.

hrfConvolve(x=grangerdata[1:100,1])
plot(hrfConvolve(grangerdata[1:100,1]))

# Compare the convolved timeseries with the raw
par(mfrow=c(2,1))
plot(x=semdata[1:100,1])
plot(hrfConvolve(x=semdata[1:100,1]))
}
\keyword{ ~hrf }
