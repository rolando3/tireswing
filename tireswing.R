require(igraph)
palette(sample(rainbow(16),16))

cMap <- function(filename,use.edge.color=FALSE,use.edge.label=FALSE,use.community=TRUE,vertex.community=walktrap.community,cex=1,font="Helvetica",title="",outputfile="") {
	df <- read.csv(filename)
	if(nrow(df)>0)
	{
		g <- graph.data.frame(df,directed=FALSE)
		E(g)$width <- log(df$cM)	   
		type="person"
		for(i in (names(df) %in% "chr")) { if(i){type="segment"} }
		if (use.edge.color) { if( type == "segment" ) { E(g)$color <- df$chr } }
		if (use.edge.label) { if( type == "segment" ) { E(g)$label <- df$chr } else { E(g)$label <- df$cM } }
		if (use.community)
		{
			com <- vertex.community(g)
			V(g)$color <- com$membership
		}
		par(mai=c(0,0,1,0))
		plot(g,vertex.label.cex=cex,edge.label.cex=cex,vertex.label.family=font,edge.label.family=font)
		title(title)
		if ( outputfile != "" )
		{
			dev.copy2pdf(file=outputfile)
		}
	}
}
mapDir <- function(path=".",numfiles=FALSE)
{
	i = 0
	for(file in list.files(path=path,pattern=".csv"))
	{
		plot.new()
		i = i + 1
		name<-substr(file,1,nchar(file)-4)
		if ( numfiles ) {
			outfile<-paste(path,"/",i,". Graph_",name,".pdf",sep="")
		} else {
			outfile<-paste(path,"/Graph_",name,".pdf",sep="")
		}
		cMap(paste(sep="",path,"/",file),title=name,outputfile=outfile)
	}
}