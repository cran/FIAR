\name{ARsem}
\alias{ARsem}
\title{Auto-regressive SEM
}
\description{Fit an auto-regressive SEM of the specified order. The function
automatically extends the given model and dataset to the model and dataset of
the specified order. The function is actually a wrapper around the function
sem() from the package lavaan.
}
\usage{
ARsem(model, data, order)
}
\arguments{
  \item{model}{A vector specifying the model of AR order 0 (AR(0)) . The vector
should be written as an n by n matrix where n is the number of regions in the
network (see example). For every expected connection ij from region i (column)
to region j (row) the vector contains '1' and '0' otherwise.
}
  \item{data}{contains all observations (rows) and variables (colums) in the
network m0. Only variables that are in the model m0 should be in the dataset.
}
  \item{order}{Integer. The order of the AR model.
}
}
\details{An AR model of order q contains the t-0, t-1, t-2,...t-q timeseries
and these timeseries are connected based on the model of order 0. The function
will transform this model in the correct AR(q) model and the data set in the
lagged data set containing all lagged variables. Let us take the very simple
example of a dataset with 2 variables X and Y. If there is an arrow from X to
Y, the function will create the AR(1) model with an additional arrow from X-1
to Y-1, from X-1 to X and from Y-1 to Y. Variables X-1 and Y-1 are
automatically created within the function.
}
\value{An object of class 'lavaanModel', for which several methods are
     available, including a 'summary' method.}
     
\references{Kim, J., Zhu, W., Chang, L., Bentler, P., and Ernst, T. (2007). Unified Structural Equation Modeling Approach for the analysis of Multisubject, Multivariate Functional MRI Data. Human Brain Mapping 85 93. 
}
\author{Bjorn Roelstraete}

\examples{
# Example dataset with three brainregions x, y, z.
head(semdata)
# Prior model with connections from (column) x to (row) y and from y to z.
model <- c(0,0,0,
           1,0,0,
	   0,1,0)
# Perform classical SEM

 fit0 <- ARsem(model,semdata)
 summary(fit0)

# Calculate AR() order of the data
ARorder(semdata,max=10)

# Compute AR(3) SEM
 fit3 <- ARsem(model,semdata,order=3)
 summary(fit3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~SEM }
\keyword{ ~auto-regressive }
