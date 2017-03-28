\name{optimiseSD}
\alias{optimiseSD}
\title{
(Spatial) optimisation of sampling designs
}
\description{
This function optimises a sampling design to achieve minimal cost for a given cost function by applying a given optimisation function. It may take into account fix locations or constrain search to a subset of locations.
}
\usage{
optimiseSD(simulations, costFun,
  locationsAll = 1:nLocations(simulations), locationsFix = integer(0),
  locationsInitial = integer(0),
  aimCost = NA, aimNumber = NA,
  optimisationFun,
  nameSave = NA, plot = FALSE, verbatim = FALSE, ...)
}

\arguments{
  \item{simulations}{ \code{\link{Simulations} object}}
  \item{costFun}{cost function, must have parameters \code{simulations} and \code{locations} (may be prepared by \code{replaceDefault} with \code{type = "costFun.optimiseSD"}); it must return cost as a single value or a list where this is the first entry}
  \item{locationsAll}{indices of the \code{locations} that are considered possible sensor locations, by default these are all locations}
  \item{locationsFix}{indices of \code{locations} with fix sensors -- to be considered when computing cost, by default this is empty}
  \item{locationsInitial}{indices of \code{locations} where sensors are initially; by default this is empty; if not indicated \code{aimNumber} random locations may be used, depending on the }
  \item{aimCost}{limit value of the \code{costFun}: the result has to fall below}
  \item{aimNumber}{if no \code{locationsInitial} given, it starts from this number of random locations}
  \item{optimisationFun}{\code{function} that executes an optimisation algorithm; must have parameters \code{simulations}, \code{costFun}, \code{locationsAll}, \code{locationsFix}, \code{locationsInitial}, \code{aimCost}, \code{aimNumber}, \code{nameSave} (use \code{type = "optimisationFun.optimiseSD"} in \code{\link{replaceDefault}}). Output must be a list with \code{SD} (vector, matrix, or list of one or several SDs) and \code{cost} (belonging to the SD); it may return other values in \code{report}}
  \item{nameSave}{\code{character} path and name (without suffix!) where to save intermediate results (sampling design of each iteration)}
  \item{plot}{if iteration is plotted (currently disabled)}
  \item{verbatim}{print intermediate results and keep sampling designs and cost of all iterations}
  \item{...}{further parameters, currently unused}
}

\value{
A list
\item{SD}{\code{list}, each entry is a matrix of the best (lowest cost) sampling designs found of a size - sampling designs in rows; including \code{locationsFix}}
\item{evaluation}{\code{data.frame}, each row belongs to the \code{SD} of one size, giving the \code{number} of sensors and the \code{cost}}
\item{aimSD}{\code{list} indicating which of the \code{SD} fulfil the given \code{aimCost} and \code{aimNumber}, can be empty if no aims given or aims not reached}
\item{report}{algorithm-specific, e.g. all tested SDs etc., see there.}
}
\details{
For examples see \code{\link{optimiseSD_genetic}}, \code{\link{optimiseSD_global}}, \code{\link{optimiseSD_greedy}}, \code{\link{optimiseSD_manual}},\code{\link{optimiseSD_ssa}}.
}
%\details{
%}
%\value{
%A list
%\item{SD}{final sampling desing}
%\item{SD_best}{best sampling design ever found}
%\item{cost}{cost of \code{SD}}
%\item{cost_best}{cost of \code{SD_best}}
%\item{report}{\code{data.frame} with columns \code{cost} and \code{accepted}, rows refer to iterations.}

%If \code{verbatim = TRUE} more information is returned:
%The \code{report} has additional columns: \code{costTest} (cost of the design tested in this iteration, no matter if it was accepted; \code{cost} refers to the currently accepted design); \code{costBest} (best cost until then); \code{chi} (random value that was used for the acceptance decision).

%In addition it returns a matrix each for \code{SDs}, \code{SDs_test}, \code{SDs_test} with the current, tested and best known sampling designs of all iterations.

%If \code{acceptanceMethod = "vanGroenigen"} it returns the cooling parameters \code{cooling} and \code{startAcceptance} as they may be computed inside the algorithm.

%It gives the iterations with \code{jumpBack} i.e. when the current sampling design was replaced by the known best.
%}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }



