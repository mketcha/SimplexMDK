require(igraph)
library(tcltk)

graph <- read.graph("~/CIS/SimplexMDK/Data/cat/mixed.species_brain_1.graphml", format = "graphml")

#tkplot(graph)

#a <- alpha.centrality(graph)

# to make directed graph undirected: as.undirected(...)

# attributes