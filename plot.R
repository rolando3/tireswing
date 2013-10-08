require(igraph)
palette(sample(rainbow(16),16))

cMap <- function(filename,type="p") {
   df <- read.csv(filename)
   if(nrow(df)>0)
   {
	   g <- graph.data.frame(df,directed=FALSE)
	   c <- walktrap.community(g)
	   E(g)$width <- log(df$cM)
	   if( type == "s") { E(g)$color <- df$chr }
	   V(g)$color <- c$membership
	   par(mai=c(0,0,1,0))
	   plot(g,vertex.label.cex=5/8,vertex.label.family="Helvetica")
   }
}

mapDir <- function(path=".",numfiles=FALSE)
{
	i = 0
	for(file in list.files(path=path,pattern=".csv"))
	{
		plot.new()
		i = i + 1
		cMap(paste(sep="",path,"/",file))
		name<-substr(file,1,nchar(file)-4)
		title(name)
		if ( numfiles ) {
			dev.copy2pdf(file=paste(path,"/",i,". Graph_",name,".pdf",sep=""))
		} else {
			dev.copy2pdf(file=paste(path,"/Graph_",name,".pdf",sep=""))
		}
	}
}