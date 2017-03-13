#####################################################
#                  Simulations-methods              #
#####################################################
# nLocations
# nPlumes
# nKinds

# modified
nLocations = function(x){
  nrow(x@values)
}
nPlumes = function(x){
  ncol(x@values)
}
nKinds = function(x){
  nlayers(x@values)
}

# old: does not work now
# nLocations = function(x){}
# nPlumes = function(x){}
# nKinds = function(x){}
# 
# nLocations.Simulations = function(x){
#   nrow(x@values)
# }
# nPlumes.Simulations = function(x){
#   ncol(x@values)
# }
# nKinds.Simulations = function(x){
#   nlayers(x@values)
# }
# 
# setMethod("nLocations", signature(x = "Simulations"), nLocations.Simulations)
# setMethod("nPlumes", signature(x = "Simulations"), nPlumes.Simulations)
# setMethod("nKinds", signature(x = "Simulations"), nKinds.Simulations)
