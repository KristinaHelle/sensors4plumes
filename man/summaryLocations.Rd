\name{summaryLocations}

\alias{summaryLocations}

\title{
Summarise values of Simulations by locations
}
\description{
This function is similar to \code{apply} for the \code{values} of \code{Simulations} location-wise. It takes a function and returns a vector with a value for each location: the result of the function applied to all values belonging to this location. This can be useful to generate maps. 
}
\usage{
summaryLocations(simulations, 
  locations = 1:nLocations(simulations), 
  plumes = 1:nPlumes(simulations), 
  kinds = 1, 
  fun, summaryFun = weightedMean, 
  weight = 1, na.rm = FALSE, ...)
}

\arguments{
  \item{simulations}{\code{Simulations}}
  \item{locations}{indices of \code{locations} to be taken into account (to apply the function to a subset of \code{simulations}); for invalid indices \code{NA} values are used}
  \item{plumes}{indices of \code{plumes} to be taken into account; for invalid indices \code{NA} values are used}
  \item{kinds}{index of kind of \code{values} to be taken into account; only one kind can be used, all but the first given index are ignored}
  \item{fun}{\code{function} to be applied; it must have a single parameter (vector) and return a single value; in addition it needs a parameter \code{na.rm} or \code{...}; possible choices are \code{sum}, \code{prod}, \code{max}, \code{min}}
  \item{summaryFun}{\code{function} to summarise the plume-wise results (as first parameter) to a global value; by default it is a weighted mean; in addition it can have the parameters \code{weight}, or use \code{...}}
  \item{weight}{to be used by \code{summaryFun}: either a numeric value of length 1 or same length as \code{plumes} or a \code{character} indicating a column of the \code{data} associated with the \code{plumes} of \code{simulations}}
  \item{na.rm}{\code{logical} how to treat missing values by \code{fun} (has no influence on \code{summaryFun})}
  \item{...}{further arguments to be forwarded to \code{fun} or \code{summaryFun}}
}

\details{
\link{summaryLocations} is similar to \code{\link{summaryPlumes}} but is not restricted to associative functions. For a more general function see \code{simulationsApply}.
To use a different \code{summaryFun} you may use \code{replaceDefault} with \code{kind = "summaryFun.summaryPlumes"} to check the function beforehand or to generate it by resetting the default parameters of an existing function.
}

\value{
A \code{list} of 
\item{summary}{the result of \code{summaryFun}}
\item{summaryLocations}{the location-wise results of \code{fun}, one value for each selected \code{location}}
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(radioactivePlumes) 
plumeSample =  c(1:5, 16:20, 30:34, 45:49)# sample of plumes: jan, apr, jul, oct

# plume detection 
# (number of plumes in plumeSample that can be detected at a level of more than 1e-7)
detection1em7 = function(x, ...){
  xExceed = x > 1e-7
  out = sum(xExceed)
  return(out)
}
plumeNr_radioactivePlumes = 
  summaryLocations(radioactivePlumes, 
                   plumes = plumeSample, fun = detection1em7, kinds = 2)
## plot map
plumeNrMap = radioactivePlumes@locations
plumeNrMap@data$plumeNr = plumeNr_radioactivePlumes[["summaryLocations"]]
spplot(plumeNrMap, zcol = "plumeNr")
}
