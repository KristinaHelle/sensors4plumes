########################################################
# optimisationCurve                                    #
########################################################

# algorithm-dependent plot of optimisation curve from results of optimiseSD
# how are fix sensors taken into account?
# postprocess
optimisationCurve = function(
  optSD, # report from optimiseSD, type depends on algorithm
  type, # ssa, genetic, greedy, complete
  nameSave, # with suffix!
  ...
  #oC,                                                 # data.frame(cost, testCost, number); rows refer to iterations
  #best,                                               # row of oC for best result
  #Imagepath,                                          # without suffix
  #TestFraction = testFraction,                        # if >0 costTest is plotted
  #kindChange = 0                                      # left of this value: grey background (SDs in grey area are without moving initial sensors)
){

  oldpar = par(no.readonly = TRUE)

  if (!missing(nameSave)){
    png(paste0(nameSave, ".png"), ...)
  }

  switch(type,
         "greedy" = {
           report = optSD$report
           # define scales
           mC = max(report[["evalSDs"]]$cost, na.rm = TRUE)
           mS = max(report[["evalSDs"]]$number, na.rm = TRUE)
           scaleRelation = mC/mS                                                         # relation of cost scale to number scale
           SscaleLag = ceiling(mS/10)                                                    # difference between two marks of nSensor scale
           Slab = SscaleLag*(0:10)                                                       # labels nSensor scale
           Clab = Slab*scaleRelation                                                     # labels cost scale

           par(mar=c(5,4,4,5) + 0.1)
           thisPlot = plot(x = 1:nrow(report[["evalSDs"]]),
                y = report[["evalSDs"]]$number,
                type = "o", col = "red",
                ylab = "", xlab = "Iteration", yaxt = "n",
                main = "Optimisation Curve", ylim = c(0, SscaleLag*11))
           usr = par("usr")
#           rect(usr[1], usr[3], kindChange + 0.5, usr[4], border = NA, col = "grey")
#            lines(x = 1:nrow(report[["evalSDs"]]),
#                  y = report[["evalSDs"]]$number,
#                  type = "o", col = "red")
           points(x = 1:nrow(report[["evalSDs"]]),
                  y = report[["evalSDs"]]$cost/scaleRelation,
                  type = "o", col = "blue")
           points(x = report$finalSDwhich,
                  y = report[["evalSDs"]]$number[report$finalSDwhich],
                  col = "red", cex = 1.5, lwd = 2, pch = 20)
           points(x = report$finalSDwhich,
                  y = report[["evalSDs"]]$cost[report$finalSDwhich]/scaleRelation,
                  col = "blue", cex = 1.5, lwd = 2, pch = 20)
# #  if(TestFraction > 0){
#     points(x = 1:nrow(report[["evalSDs"]]),
#            y = report[["evalSDs"]]$testCost/scaleRelation,
#            type = "o", col = "cyan")
#     points(x = best,
#            y = report[["evalSDs"]]$testCost[report$finalSDwhich]/scaleRelation,
#            col = "cyan", cex = 1.5, lwd = 2, pch = 16)
#     mtext("(test)", 2, col = "cyan", las = 1, line = 1.5, adj = 0,
#           at = SscaleLag*10.5, cex = 0.6)
#   }
           axis(4, Slab, labels = Slab,
                col.axis = "red",
                las = 1, cex.axis = 0.6, adj = 0,
                line = 0, mgp = c(0.1, 0.3, 0), tck = -0.01)
           mtext("Number of Sensors",
                 4, col = "red", line = 2)#, adj = 1, at = SscaleLag*11, cex = 0.6, tck = -0.01)
           axis(2, Slab,
                labels = signif(Slab*scaleRelation, digits = 3),
                col.axis = "blue", las = 1, cex.axis = 0.6, mgp = c(0.1, 0.3, 0),tck = -0.01)
           mtext("Cost", 2, col = "blue", line = 2)#, adj = 0, at = SscaleLag*11, cex = 0.6)
         },
         "ssa" = {
           report = optSD$report
           if (nrow(report$iterations) > 1){
             report$iterations = report$iterations[1:(nrow(report$iterations) -1),]  
           }
           
           if (all(is.element(c("costTest", "cost", "costBest"), names(report$iterations)))){
             thisPlot = plot(report$iterations$costTest,
                             ylim = range(report[,c("costTest", "cost", "costBest")]),
                             pch = ".", col = 2,
                             xlab = "Iteration", ylab = "Cost", main = "Optimisation Curve")
             lines(report$iterations$cost, col = 4)
             lines(report$iterations$costBest, col = 3)
           } else {
             thisPlot = plot(report$iterations$cost,
                             pch = ".",
                             xlab = "Iteration", ylab = "Cost", main = "Optimisation Curve")
           }
         },
         "genetic" = {
           report = optSD$report
           par(mfrow = c(1,2))
           # optimisation curve
           thisPlot = plot(report$best,
                           ylim = range(report$best, report$mean),
                           col = 3,
                           xlab = "Iteration", ylab = "Cost", main = "Optimisation Curve")
           lines(report$mean, col = 4)
           # population
           nSensors = apply(FUN = sum, X = report$population, MARGIN = 1)
           thisPopulation = plot(x = nSensors, y = report$evaluations,
                                 ylab = "Cost", xlab = "Number of Sensors", main = "Population")
         },
         "global" = {
           # which sensor detects how many plumes?
           optSD
           detected = matrix(nrow = nrow(optSD$SD[[1]]), ncol = ncol(optSD$SD[[1]]))
           for (i in 1:nrow(optSD$SD[[1]])){
             for (j in 1:ncol(optSD$SD[[1]])){
               detected[i,j] =
                 sum(optSD$report$detectable[optSD$SD[[1]][i,1:j],])
               if (j > 1){
                 detected[i,j] = detected[i,j] - sum(detected[i, 1:(j - 1)])
               }
             }
           }
           barplot(t(detected), xlab = "Sampling Designs", ylab = "Detected Plumes")
         },
         "manual" = {
           report = optSD$report
           nSensors = sapply(FUN = length, X = report$locationsCurrent)
           # define scales
           mC = max(report$cost[1:length(report$locationsCurrent)], na.rm = TRUE)
           mS = max(nSensors, na.rm = TRUE)
           scaleRelation = mC/mS                                                         # relation of cost scale to number scale
           SscaleLag = ceiling(mS/10)                                                    # difference between two marks of nSensor scale
           Slab = SscaleLag*(0:10)                                                       # labels nSensor scale
           Clab = Slab*scaleRelation                                                     # labels cost scale

           par(mar=c(5,4,4,5) + 0.1)
           thisPlot = plot(x = 1:length(report$locationsCurrent),
                           y = nSensors,
                           type = "o", col = "red",
                           ylab = "", xlab = "Iteration", yaxt = "n",
                           main = "Optimisation Curve", ylim = c(0, SscaleLag*11))
           usr = par("usr")
           points(x = 1:length(report$locationsCurrent),
                  y = report$cost[1:length(report$locationsCurrent)]/scaleRelation,
                  type = "o", col = "blue")
           axis(4, Slab, labels = Slab,
                col.axis = "red",
                las = 1, cex.axis = 0.6, adj = 0,
                line = 0, mgp = c(0.1, 0.3, 0), tck = -0.01)
           mtext("Number of Sensors",
                 4, col = "red", line = 2)#, adj = 1, at = SscaleLag*11, cex = 0.6, tck = -0.01)
           axis(2, Slab,
                labels = signif(Slab*scaleRelation, digits = 3),
                col.axis = "blue", las = 1, cex.axis = 0.6, mgp = c(0.1, 0.3, 0),tck = -0.01)
           mtext("Cost", 2, col = "blue", line = 2)#, adj = 0, at = SscaleLag*11, cex = 0.6)
         })
  par(oldpar)
  if (!missing(nameSave)){
    dev.off()
  }
}

