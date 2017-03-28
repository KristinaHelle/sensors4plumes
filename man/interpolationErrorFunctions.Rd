\name{interpolationErrorFunctions}
\alias{absError}
\alias{delineationError}
\alias{absErrorMap}
\alias{delineationErrorMap}
\alias{interpolationErrorFunctions}

\title{
Compare original and interpolated maps point-wise
}
\description{
\code{absError} and \code{delineationError} compare an interpolated to an original map (the values in each point respectively). They are to be used as \code{fun_error} in \code{interpolationError}. \code{absError} is the absolute difference between the values. \code{delineationError} compares areas above a given threshold, and indicates false positive and false negative classification.

\code{absErrorMap} and \code{delineationErrorMap} are to be used as \code{fun_l} in \code{interpolationError}, they compute the average over all plumes of the respective error functions to generate a common map.
}
\usage{
absError(x, nout = 1)
delineationError(x, nout = 1, threshold = 1e-7, weightFalseNeg = 5)
absErrorMap(x, nout = 1)
delineationErrorMap(x, nout = 1, threshold = 1e-7, weightFalseNeg = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of length 2 with \code{x[1]} the value and \code{x[2]} the interpolated value }
  \item{nout}{length of output -- needed as the function is to be called via \code{simulationsApply} as \code{fun_pl} or \code{fun_l}}
  \item{threshold}{threshold to classify all locations}
  \item{weightFalseNeg}{weight for false negative classification, false positive classification is weighted 1}
}

\value{
\code{absError} returns a single numeric value, the absolute difference. 

\code{delineationError} first computes if the original is above the threshold and if the interpolated is above the threshold; then it determines the result: 0 if classifications agree, 1 for false positive (i.e. original is below threshold, interpolated is above), and \code{weightFalseNeg} for false negative.
%returns three values: if the original is above the threshold; if the interpolated is above the threshold; the result: 0 if classifications agree, 1 for false positive (i.e. original is below threshold, interpolated is above), and \code{weightFalseNeg} for false negative.

\code{absErrorMap} and \code{delineationErrorMap} return single values.
}
%\references{
%}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
## see interpolationError
}
