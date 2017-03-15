\name{optimisationCurve}
\alias{optimisationCurve}

\title{Plot Optimisation Curve}

\description{
Generates algorithm-specific plots of optimisation or optimal sampling design(s).
}
\usage{
optimisationCurve(report, type, nameSave, ...)
}
\arguments{
  \item{report}{ \code{report} from a result of \code{optimiseSD}, algorithm-specific format}
  \item{type}{character, indicating type of report, of \code{"ssa", "genetic", "greedy", "global", "manual"}}
  \item{nameSave}{path of file where to save, without suffix, generates a .png file}
  \item{...}{parameters to be forwarded to \code{\link[grDevices]{png}}}
}
\details{
The type of plot depends on the algorithm: 

\code{"ssa"}: curve of the cost in each iteration: proposed design (red dot), accepted design (blue line), best design until now (green line)

\code{"genetic"}: two plots in one panel. Optimisation curve: mean (blue line) and best (green circles) cost in each iteration. Population: cost versus number of sensors in the final population. It may be useful to choose larger \code{width} as plots are placed beside each other. 

\code{"greedy"}: cost (blue) and number of sensors (red) in each iteration, combined in one plot with adjusted scales. Lowest cost is marked by solid dots.

\code{"global"}: a barplot of all sampling designs, showing how many plumes are detected by 1st, 2nd, etc. sensor. 

\code{"manual"}: similar to the plots of \code{"greedy"}.
}
\value{
Generates a plot or a file with it, no value returned.
}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de} }

\examples{
data(SDgreedy)
curve_greedy1 = optimisationCurve(
  report = SDgreedy,
  type = "greedy")

data(SDgenetic)
curve_genetic1 = optimisationCurve(
  report = SDgenetic$report,
  type = "genetic")

data(SDglobal)
curve_global1 = optimisationCurve(
  report = SDglobal,
  type = "global")
  
data(SDmanual)  
curve_manual1 = optimisationCurve(
  report = SDmanual,
  type = "manual")

\dontshow{
\dontrun{
# error
data(SDssa)
curve_ssa1 = optimisationCurve(
  report = SDssa$report,
  type = "ssa")
}
}
  
\dontrun{
# generate a file
curve_global1 = optimisationCurve(
  report = SDglobal,
  type = "global",
  nameSave = "optSD_global",
  width = 600, height = 300)
}  

}
