\name{DCMex}
\alias{DCMex}
\docType{data}
\title{
Demo dataset
}
\description{
A toy dataset containing all parameters of an estimated DCM
}
\usage{data(DCMex)}
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
    
    \item{\code{HPF}}{Length of High Pass filter in seconds}
    \item{\code{x}}{Number of states (5 per region)}
    \item{\code{sf}}{Stimulus function of the experiment. One column per experimental condition}
    \item{\code{s}}{Information about stimulusfunction}
    \item{\code{T0}}{Information about stimulusfunction}
    \item{\code{dt0}}{Information about stimulusfunction}
    
    \item{\code{X0}}{Confounds or null space}
    \item{\code{sim}}{Simulated timeseries from model with SNR=1 and ar=0}
    \item{\code{priors}}{List containing prior parameter covariances (pC), parameter expectations (pE), hemodynamic expectations (qE), and hemodynamic covariances (qC)}
    \item{\code{Ep}}{Posterior model parameters}
    
    \item{\code{A}}{Posterior anatomical connections}
    \item{\code{B}}{Posterior functional connections}
    \item{\code{C}}{Posterior input connections}
    \item{\code{H}}{Posterior hemodynamic parameters}
    \item{\code{Cp}}{Posterior parameter covariances}
    \item{\code{F}}{Log evidence}
    \item{\code{Ce}}{Posterior error covariances}
    \item{\code{names}}{Names of regions} 
}
}
\examples{
names(DCMex)
}

