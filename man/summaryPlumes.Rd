\name{summaryPlumes}

\alias{summaryPlumes}
\alias{weightedMean}

\title{
Summarise values of Simulations by plumes
}
\description{
This function is similar to \code{apply} for the \code{values} of \code{Simulations} plume-wise. It takes a function and returns a vector with a value for each plume: the result of the function applied to all values belonging to this plume. As \code{values} may be too big to keep in memory, the function is applied chunk-wise; therefore it only works for associative functions. A (slower) method without this restriction is %\code{\link{simulationsApply}}.
It is combined with a global summary function and thus provides a general basic form of a cost function.
}
\usage{
summaryPlumes(simulations,
  locations = 1:nLocations(simulations), 
  plumes = 1:nPlumes(simulations), 
  kinds = 1, 
  fun, summaryFun = weightedMean,
  weight = 1, na.rm = FALSE, ...)
weightedMean(x, weight, na.rm = FALSE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simulations}{\code{Simulations}}
  \item{locations}{indices of \code{locations} to be taken into account (to apply the function to a subset of \code{simulations}); for invalid indices \code{NA} values are used; multiple are ignored}
  \item{plumes}{indices of \code{plumes} to be taken into account; for invalid indices \code{NA} values are used}
  \item{kinds}{index of kind of \code{values} to be taken into account; only one kind can be used, all but the first given index are ignored}
  \item{fun}{\code{function} to be applied; it must have a single parameter (vector), be associative, i.e. \code{fun(c(a,b)) = fun(c(fun(a),b))}, work also on vectors of length 0 and return a single value; in addition it needs a parameter \code{na.rm} or \code{...}; possible choices are \code{sum}, \code{prod}, \code{max}, \code{min}}
  %\item{initial}{initial value for the \code{fun}, should be neutral i.e. \code{fun(a) = fun(c(initial, a))} for all values; e.g. \code{Inf} for \code{fun = "min"} (if \code{fun} is a \code{character}, the appropriate \code{initial} is used automatically (it may be overwritten by the user-defined value)}
  \item{summaryFun}{\code{function} to summarise the plume-wise results (as first parameter) to a global value; by default it is a weighted mean; in addition it can have the parameters \code{weight}, or use \code{...}}
  \item{weight}{to be used by \code{summaryFun}: either a numeric value of length 1 or same length as \code{plumes} or a \code{character} indicating a column of the \code{data} associated with the \code{plumes} of \code{simulations}}
  \item{na.rm}{\code{logical} how to treat missing values by \code{fun} (has no influence on \code{summaryFun})}
  \item{...}{further arguments to be forwarded to \code{fun} or \code{summaryFun}}
  \item{x}{vector}
}

\details{
To use a different \code{summaryFun} you may use \code{replaceDefault} with \cr \code{kind = "summaryFun.summaryPlumes"} to check the function beforehand or to generate it by resetting the default parameters of an existing function.
}

\value{
A \code{list} of 
\item{summary}{the result of \code{summaryFun}}
\item{summaryPlumes}{the plume-wise results of \code{fun}, one value for each selected \code{plume}}
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(radioactivePlumes) 
# sample of locations, e.g. proposed sensor set
locSample = sample.int(nLocations(radioactivePlumes), 10) 

# plume detection 
# (number of plumes not detected at any of locSample at a level of more than 1e-7)
plumeMin_radioactivePlumes = 
  summaryPlumes(radioactivePlumes, locations = locSample, kinds = 2, 
  fun = "min",
  summaryFun = function(x, weight = 1, na.rm = TRUE){sum(x < 1e-7)})
}
