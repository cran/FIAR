\name{dcmParam}
\alias{dcmParam}

\title{DCM object builder}
\description{Automated step by step procedure to enter parameters needed for
the DCM analysis.  Everything is stored in an DCM list that can be used to
generate timeseries or estimate the model.  The function also immediately
calculates the High pass filter (HPF) and stimulus (SF) function. 
}
\usage{
dcmParam(a, b, c, ons = list(), dur = list(), v, n, m, TR, 
          h = c(0.65, 0.41, 0.98, 0.32, 0.34, 0), names=c(),
          TE = 0.04, T = 16, x = 5 * n, HPF=0, auto = FALSE)
}

\arguments{
\item{a}{vector of length n*n representing the anatomical connections between
brain regions. The 'a' vector should look like an n*n matrix with a non zero
element (in Hz.) on location ij when there is an expected connection from
region j (column) to region i (row)(See examples).}
\item{b}{vector of length n*n*m representing the functional connections. Should
be written as m n*n matrices whith non zero elements ijk (in Hz.) if input k
influences the connection from region j (column) to region i (row) and zero
otherwise.}
\item{c}{vector of length n*m representing the input connections. Should be
written as an m by n matrix whith non zero elements kj (in Hz.) if input k
(row) influences region j (columns).}
\item{ons}{list containing vector of onsets (in scans) for every input.}
\item{dur}{list containing vector of durations (in scans) for every input.
dur=0 represents event-related inputs. If durations are equal over the entire
experiment, the number only needs to be entered once. If durations differ, the
length of the durationvector should match the length of the onsetvector.}
\item{v}{number of scans}
\item{n}{number of regions in model}
\item{m}{number of inputs (experimental conditions)}
\item{TR}{repetition time of the experiment}
\item{h}{ the parameters of the hemodynamic model
 % \item{H[1]}{signal decay} \item{H[2]}{autoregulation} \item{H[3]}{transit
 % time} \item{H[4]}{exponent for Fout(v)} \item{H[5]}{resting oxygen
 % extraction} \item{H[6]}{ratio of intra- to extra-vascular components of the
 % gradient echo signal                   %}
  }
\item{names}{names of brain regions}
\item{TE}{echo time in seconds (default = 0.04 s)}
\item{T}{number of timebins (default = 16)}

\item{x}{number of states (= 5 times number of regions)}
\item{HPF}{High pass filter in seconds (default = 0 seconds)}
\item{auto}{logical. If FALSE (default) the prespecified DCM object is used and
the function only serves to construct the HPF and SF. If TRUE, the DCM object
need not be prespecified and is constructed in a step by step procedure,
whereafter the HPF and SF are constructed. 
}
}
\value{DCM list containing all above mentioned model and scanner parameters.
}
\author{Bjorn Roelstraete}
\examples{
## Specify connectivity parameters in a 3 region network with connections 
## from region 1 to region 2 with a strenght of .8 Hz and region 2 to 
## region 3 with .65Hz.
a <- c( 0, 0, 0,
       .7, 0, 0,
        0,.4, 0)   

## Specify 2 experimental manipulations (inputs) where the first directly
## influences region 1 with .4 Hz. and the second region 2 with .2 Hz
c <- c(.4, 0, 0,
        0,.5, 0) 

## Specify the functional connectivities between region 1 and 3 of .2 Hz. 
## caused by input 1. Input 2 influences the functional connectivity from 
## region 3 to region 2.
b <- c(0, 0, 0,
       0, 0, 0,
       0,.2, 0,
       0, 0, 0,
       0, 0,.2,
       0, 0, 0)		

v <- 240            # number of scans 
n <- 3              # number of regions
m <- 2              # number of inputs = number of rows in DCM$c

## The onsets of input 1 are at scan 0, 30, 60, 120, and 200. The onsets of
## input 2 at scan 30, and 120.

ons.input1 <- c(0, 60, 120, 180)
ons.input2 <- c(0, 30, 60, 90, 120)

## The duration of input 1 is always 30 scans. The duration of input 2 is 15
## scans.

dur.input1 <- 30
dur.input2 <- 15
TR <- 1
DCM <- dcmParam(a, b, c, ons=list(ons.input1,ons.input2), 
                 dur=list(dur.input1,dur.input2), v, n, m, TR=1)
names(DCM)

# Or construct DCM by automated step by step:
# DCM <- dcmParam(auto=TRUE)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~DCM }
