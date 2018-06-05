\name{grangerdata}
\alias{grangerdata}
\docType{data}
\title{
Demo dataset
}
\description{
A toy dataset containing 5 autoregressive timeseries generated from the model by Baccala and Sameshima (Biol. Cybern. 2001).
}
\usage{data(grangerdata)}
\format{
A data frame of 2000 observations of 3 variables.
\describe{
    \item{\code{x}}{Time series at region 1}
    \item{\code{y}}{Time series at region 2}
    \item{\code{z}}{Time series at region 3}
  \item{\code{q}}{Time series at region 4}
  \item{\code{w}}{Time series at region 5}
  }
}
\examples{
head(grangerdata)
}

