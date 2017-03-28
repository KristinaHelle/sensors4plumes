\name{replaceDefault}
\alias{replaceDefault}
\title{
Replace default values in functions and check parameter list}
\description{
\code{replaceDefault} takes a function and a list of values to replace the default values. In addition it checks that the resulting function is fit for use in other functions, e.g. as cost function in \code{optimiseSD}. }
\usage{
replaceDefault(fun, newDefaults = list(), type)
}
\arguments{
  \item{fun}{
 function }
  \item{newDefaults}{
 \code{\link{list}} of new default values: elements must be named, names have to be those of parameters of \code{fun} (not necessarily all of them); then default values of \code{fun} are replaced; elements without or with other names cause a warning. }
  \item{type}{
  character string indicating the purpose of the resulting function, e.g. to be used as \code{summaryFun} in \code{summaryPlumes}; these function calls require certain parameters: it is checked if the \code{fun} has them and if all additional parameters have default values (after replacing defaults), and generates warnings, if not. For possible values, see details.}
}

\details{
Dependent on the \code{type}, a list of required parameters is generated and compared to the parameters of \code{fun}: all required parameters must exist and all other parameters need to have default values (in the resulting function, when defaults of \code{newDefaults} have been inserted)

"summaryFun.summaryPlumes" (\code{summaryFun} in \code{summaryPlumes}): "x", "weight" \cr
"fun.simulationsApply" (\code{fun}, \code{fun_l}, \code{fun_p}, \code{fun_pl} in \code{simulationsApply}): "x", "nout" \cr
"funR.simulationsApply" (\code{fun_Rp} or \code{fun_Rl} in \code{simulationsApply}): "x", "nout", "weight" \cr
"funRR.simulationsApply" (\code{fun_Rpl} in \code{simulationsApply}): "x", "nout", "weight_l", "weight_p" \cr
"fun_interpolation.interpolate" (\code{fun_interpolation} in \code{interpolate}): "data", "newdata", "y" \cr
%"fun_interpolationSplit.interpolate" (\code{fun_interpolation} in \code{interpolate} with extra parameters for regional interpolation): "data", "newdata", "y", "dataLoc", "newdataLoc" \cr
"fun.spatialSpread" (\code{fun} in \code{spatialSpread}): "allLocations", "locations"\cr
"fun_R.spatialSpread" (\code{fun_R} in \code{spatialSpread}): "x" \cr
"costFun" (\code{costFun} in \code{optimiseSD}): "simulations", "locations" \cr
"optimisationFun.optimiseSD" (\code{optimisationFun} in \code{optimiseSD}): "simulations", "costFun", "locationsAll", "locationsFix", "locationsInitial", "aimcost", "aimNumber", "nameSave" \cr
"evalFunc.rbga.bin"(\code{evalFunc} in \code{rbga.bin}): "chromosome"
%"rbga.bin" (?) "size", "suggestions", "popSize", "iters", "mutationChance", "elitism", "zeroToOneRatio", "monitorFunc", "evalFunc", "showSettings", "verbose"

For \code{"fun.simulationsApply"} and \code{"funR.simulationsApply"} the parameter "nout" must exist and have   default value.

%\code{optimiseSD} takes the parameter \code{cost}; this has to be a function with a parameter \code{locations}, and maybe additional parameters, that must have default values. Similarily it uses functions for the optimisation algorithm (e.g. \code{spatialSimulatedAnnealing}), for these functions all but some parameters must have default values. 
\code{replaceDefault} checks, if \code{fun} is a function. It does not check if the new defaults have the same class as the old ones. 
}
\value{
A list
\item{fun}{function \code{fun} with replaced default values}
\item{accept}{\code{logical} if the function passed all tests (i.e. no warnings) -- mainly needed when \code{replaceDefault} is called inside of other functions (e.g. \code{simulationsApply})}
}

%\references{
%% ~put references to the literature/web site here ~
%}
\author{ Kristina B. Helle, \email{kristina.helle@uni-muenster.de}}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# examples with warnings
fun1 = "a"
fun_1 = replaceDefault(fun1)
fun2 = function(x){sum(x)}
fun_2 = replaceDefault(fun2, type = "fun.simulationsApply")
fun3 = function(simulations, locations){1}
# change default of 'simulations' 
fun_3 = replaceDefault(fun3, type = "costFun", 
  newDefaults = list(nout = 1, simulations = 10))  
# example without warnings 
# (new default values have different class, not tested)
fun6 = function(x = 17, weight = "a", extraWeight = matrix(1:12, nrow = 3)){}
fun_6 = replaceDefault(fun6, 
  newDefaults = list(weight = 3, extraWeight = "c"), 
  type =  "summaryFun.summaryPlumes") 
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
