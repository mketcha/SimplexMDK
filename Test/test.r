require(igraph)
library(tcltk)
#library(FlashGraphR)

## igraph
graph <- read.graph("~/Simplex/SimplexMDK/Data/Test/cat/mixed.species_brain_1.graphml", format = "graphml")
graph <- as.undirected(graph)
a <- alpha.centrality(graph)

A = adjacency.spectral.embedding(graph,20)
#tkplot(graph)



# to make directed graph undirected: as.undirected(...)

# attributes

## FlashGraphR
#fg2 <- fg.load.igraph(graph)
#fg2 <- fg.load.graph("~/Simplex/SimplexMDK/Data/Test/wikiVote.txt")
#res <- fg.clusters(fg, "weak")