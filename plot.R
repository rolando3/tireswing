require(igraph)
palette(sample(rainbow(16),16))

cMap <- function(filename) {
   df <- read.csv(filename)
   g <- graph.data.frame(df,directed=FALSE)
   c <- walktrap.community(g)
   E(g)$width <- log(df$cM)
   V(g)$color <- c$membership
   par(mai=c(0,0,1,0))
   plot(g,vertex.label.cex=5/8,vertex.label.family="Helvetica")
}

mapDir <- function(path=".")
{
	i = 0
	for(file in list.files(path=path,pattern=".csv"))
	{ 
		i = i + 1
		cMap(file)
		name<-substr(file,1,nchar(file)-4)
		title(name)
		quartz.save(type="pdf",paste(path,"/",i,". Graph_",name,".pdf",sep=""))
	}
}