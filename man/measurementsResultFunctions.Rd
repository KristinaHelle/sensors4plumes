\name{mesurementsResultFunctions}

\alias{measurementsResultFunctions}
\alias{singleDetection}
\alias{multipleDetection}
\alias{earlyDetection}

\title{
Cost functions by plume-wise summary of values at locations 
}
\description{
\code{singleDetection} is average total dose of non-detected plumes. 

\code{multipleDetection} is the average 
fraction of possible additional detections: given a plume can be detected in 10 locations and the given sensors detect it in 3, then the cost would be 7/10; this is averaged over all plumes.

\code{earlyDetection} is the average time between the plume entering the area and its detection.
}

\usage{
singleDetection(simulations, locations, plot = FALSE)
multipleDetection(simulations, locations, plot = FALSE)
earlyDetection(simulations, locations, plot = FALSE)  
}
\arguments{
  \item{simulations}{\code{Simulations} object, must contain the required data, see details.}
  \item{locations}{indices of \code{locations}, sensor positions}% multiple, invalid, order?
  \item{plot}{if a map is generated (takes time), only implemented yet for \code{singleDetection}, else it is automatically set to \code{FALSE}}
}
%\details{
%For these cost functions assume that for a set of plumes \eqn{p \in P} the dose rate \eqn{r_p(x,t)} is known for all locations of interest \eqn{x \in X} and time steps \eqn{t \in T}. A plume is detected when \eqn{r_p(x,t) \geq r_{threshold}} for a threshold dose rate \eqn{r_{threshold}}. For the first two functions we only take into account if plumes are detected (exceed a threshold \eqn{r_{threshold}}) at all, no matter when: 
%\deqn{
%I_x(p) = \left{
%  \begin{array}{ll}
%  1 & \mbox{if $\max_{t \in T}\{r_p(x,t)\} \geq r_{threshold}$} \\
%  0 & \mbox{else ($p$ never detected at $x$).} \\
%  \end{array} \right.}
%}

%\code{singleDetection} is average total dose of non-detected plumes
%\deqn{
%  c_{sd}(S) = 1/nP \sum_{p \in P} (1 - I_S(p)) \cdot D(p) 
%}
%where \eqn{P} are all plumes (\eqn{nP} their number), \eqn{D(p)} is the total dose \eqn{D(p) = \sum_{x \in X}\sum_{t \in T} r_p(x,t)}
%and \eqn{I_S(p) = 1 - \prod_{s \in S}(1 - I_s(p))} indicates if the plume can be detected in any of the sensor locations \eqn{s \in S \subset X}. 


%\code{multipleDetection} only takes into account in how many locations a plume could be detected and how many sensors detect it: 
%\deqn{
%c_{md}(S) = 1/nP \sum_{p \in P}  (1 - \fraq{\sum_{s \in S} I_s(p)}{\sum_{x \in X} I_x(p)})
%}

%\code{earlyDetection} is the average time until detection (compared to the first possible detection in this area). It is based on the first detection in each location 
%\deqn{
%t0_x(p) = \min\{t \in T: r_p(x,t) \geq r_{threshold}\}.
%}
%Then the first time a plume can be detected anywhere is the minimum \eqn{t0(p) = \min_{x \in X} \{t0_x(p)\}} and the cost is the time difference between $t0(p)$ and earliest detection at any of the sensors
%\deqn{
%c_{ed}(S) = 1/nP \sum_{p \in P} 1/t_{max}(t0(p) - \min_{s \in S} \{t0_s(p)\}).
%}
%For calibration it is divided by $t_{max}$ the time assumed until detection in cases where plumes are not detected at all; we use twice the longest occuring time until detection \eqn{t_max = 2 * \max_{x \in X}\{t0_x(p)\} = 2 * 6.0480e+05}.


%The functions are tailored to the kind of data in \code{radioactivePlumes_area} and \code{radioactivePlumes_local} as the \code{values} there are: \cr
%\code{"maxdose"}: \eqn{\max_{t \in T}r_p(x,t)} (maximal dose during plume passage, can be used to compute detection \eqn{I_x(p)})\cr
%\code{"finaldose"}: \eqn{d_x(p) = \sum_{t \in T} r_p(x,t)} (cummulative dose, can be used to compute total dose per plume \eqn{D(p)})\cr
%\code{"time"}: \eqn{t0_x(p)} (time until first detection at a threshold of \eqn{t_{threshold} = 1e-7 Sv/h})

%The three cost functions can be computed from these values. As we assume that cost will be computed several times for different sensors (as in optimisation) some intermediate results have to be computed beforehand, the code is given in the examples:
%for \code{singleDetection} and \code{multipleDetection} the \code{values} need a layer \code{detecable} (\eqn{I_x(p)}) \cr
%for \code{singleDetection} the \code{plumes} needs a column \code{totalDose} (\eqn{D(p)}) with the total dose of the plume in the area of interest as weight \cr
%for \code{multipleDetection} the \code{plumes} needs a column \code{nDetectable} (\eqn{\sum_{x \in X} I_x(p)}) with the number of locations where this plume can be detected \cr
%for for \code{earlyDetection} the \code{plumes} needs a column \code{earliestDetection} (\eqn{t0(p)}) which is the time until first detection in any location.
%}

\value{
  \code{"cost"}: global cost
  \code{"costPlumes"}: plume-wise intermediate result, not to be interpreted as plume-specific cost
}

\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
# load data
demo(radioactivePlumes_addProperties)

# define 'min' function without warning for empty sets
min_ = function(x, na.rm = TRUE){
  if (length(x) == 0){ 
    out = Inf
  } else {
    out = min(x, na.rm = na.rm)
  } 
  return(out)
} 

# preprocess data to provide required input  
### earliest possible detection of plume inside the area
radioactivePlumes@plumes$earliestDetection = 
  summaryPlumes(radioactivePlumes, fun = min_, kinds = "time", na.rm = TRUE)[[2]]
  
# sample locations (sensors)
sampleLocations1 = sample.int(nLocations(radioactivePlumes), 10)

# compute cost
singleDetection1 = singleDetection(
  simulations = radioactivePlumes,
  locations = sampleLocations1)
  
multipleDetection1 = multipleDetection(
  simulations = radioactivePlumes,
  locations = sampleLocations1)
  
earlyDetection1 = earlyDetection(
  simulations = radioactivePlumes,
  locations = sampleLocations1)
}
