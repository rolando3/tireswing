require(igraph)
palette(sample(rainbow(16),16))

cMap <- function(filename) {
   df <- read.csv(filename)
   g <- graph.data.frame(df,directed=FALSE)
   c <- walktrap.community(g)
   E(g)$width <- log(df$cM)
   V(g)$color <- c$membership
   plot(g,vertex.label.cex=5/8,vertex.label.family="Helvetica")
}

