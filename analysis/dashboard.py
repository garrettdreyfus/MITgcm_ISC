from graph import folderMap 
import graph
import matplotlib.pyplot as plt


runsdict = {\
        "icefront-64":{"specialstring":['y100','y250'], "marker":["^","^"] ,"color":["orange","purple"],"description":["Different ice front distances"]},\
        "shelfdepth-16":{"specialstring":['d500','d600','d700'], "marker":["x","x","x","x"] ,"color":["purple","gray","green","red"],"description":["Different shelf depths"]},\
        "widthexp-GLIB-explore-32":{"specialstring":['w50','w100','w250'], "marker":["o","o","o"] ,"color":["orange","purple","gray"],"description":["Different shelf widths" ]},\
        "slope-22":{"specialstring":["s300","s150","s0"], "marker":["D","D","D"] ,"color":["red","green","blue"],"description":["more steep slope"]},\
        "reference":{"specialstring":["ref"], "marker":["p"] ,"color":["black"],"description":["more steep slope"]}\
        }

#generateRunsTable(runsdict)
#legendFunction(runsdict)
#crossSectionAnim("/home/garrett/Projects/MITgcm_ISC/experiments/widthexp-GLIB-explore-32/at0w250/results/","")

#fastExplosionCheck(runsdict)
#folderMapRefresh(runsdict)
#folderMapTimeSeries(runsdict)
#folderMap(runsdict)
#plt.show()
#folderMapGeneric(gprimeWidth,runsdict)
#plt.show()
#plt.show()
#crossSectionAverage("/home/garrett/Projects/MITgcm_ISC/experiments/reference/at125/results","Reference")
graph.crossAndMelt("/home/garrett/Projects/MITgcm_ISC/experiments/reference/at125/results","Reference")
plt.show()
#graph.circulationFigure("/home/garrett/Projects/MITgcm_ISC/experiments/reference/at125/results","Reference")
#plt.show()

#folderMapGeneric(steadyStateAverageSimple,runsdict)
#folderMapGeneric(gprimeWidth,runsdict)
#folderMapGeneric(saltBudget,runsdict)
#plt.show()
#folderMapGeneric(steadyStateAverageSimple,runsdict)
#plt.show()

