###############################################################
# test postprocessing                                         #
###############################################################
# use results form optimiseSD of all different algorithms

# report = optGreedy5
curve_greedy1 = optimisationCurve(
  report = optGreedy5,
  type = "greedy")
curve_greedy2 = optimisationCurve(
  report = optGreedy5,
  type = "greedy",
  nameSave = "curve_greedy1.png",
  width = 300, height = 400, bg = "yellow")

# report = optSSA1$report
curve_ssa1 = optimisationCurve(
  report = optSSA1$report,
  type = "ssa")
curve_ssa2 = optimisationCurve(
  report = optSSA1$report,
  type = "ssa",
  nameSave = "curve_ssa1.png")

# report = optSD_gen1$report
curve_genetic1 = optimisationCurve(
  report = optSD_gen1$report,
  type = "genetic")
curve_genetic2 = optimisationCurve(
  report = optSD_gen1$report,
  type = "genetic",
  nameSave = "curve_genetic1.png")

# report = optSD_global2
curve_global1 = optimisationCurve(
  report = optSD_global2,
  type = "global")
curve_global2 = optimisationCurve(
  report = optSD_global2,
  type = "global",
  nameSave = "curve_global1.png")

# report = optSD_manual2
curve_manual1 = optimisationCurve(
  report = optSD_manual2,
  type = "manual")
curve_manual2 = optimisationCurve(
  report = optSD_manual2,
  type = "manual",
  nameSave = "curve_manual1.png")
