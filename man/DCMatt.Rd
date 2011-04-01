\name{DCMatt}
\alias{DCMatt}
\docType{data}
\title{
Demo dataset
}
\description{
An DCM object containing all parameters for estimating the attention to visual motion DCM from the SPM8 manual.
}
\usage{data(DCMatt)}
\format{
A list containing:
\describe{
    \item{\code{a}}{Prior anatomical connections}
    \item{\code{b}}{Prior functional connections}
    \item{\code{c}}{Prior input connections}
    \item{\code{h}}{Prior hemodynamic parameters}
    \item{\code{ons}}{List containing onsets per experimental conditions in scans}
    \item{\code{dur}}{List containing durations per experimental conditions in scans}
    \item{\code{T}}{Number of timebins}
    \item{\code{TR}}{repetition time in seconds}
    \item{\code{TE}}{Echo time in seconds}
    \item{\code{m}}{Number of inputs}
    \item{\code{v}}{Number of scans}
    \item{\code{n}}{Number of regions}
    \item{\code{names}}{names of regions}
    \item{\code{x}}{Number of states (5 per region)}
    \item{\code{X0}}{Confounds or null space}
    \item{\code{names}}{Names of regions} 
}
}
\examples{
names(DCMatt)
}

