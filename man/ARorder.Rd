\name{ARorder}
\alias{ARorder}
\title{Estimate AR order of Multivariate timeseries
}
\description{Compute AIC of a series of Multivariate AR models and returns the order of the model which minimizes this AIC.
}
\usage{
ARorder(data, min=1,max = 10,type='AIC')
}
\arguments{
  \item{data}{timeseries to compute autoregressive order of.
}
\item{min}{Minimum order of AR model to check.
}  
\item{max}{Maximum order of AR model to check.
}
\item{type}{Use AIC or BIC to compute model order.
}
}
\value{returns the order (>=1) of the autoregressive model which minimizes the AIC or BIC
}
\author{Bjorn Roelstraete}

\examples{
# Compute the AR order of semdata based on AIC, with a maximum order of 10 to reduce computing time.
ARorder(semdata,max=10)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~timeseries }
