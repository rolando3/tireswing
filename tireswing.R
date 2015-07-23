require(igraph)

cMap <- function(filename,chromosome=0,use.edge.color=FALSE,use.edge.label=FALSE,use.community=FALSE,vertex.community=walktrap.community,cex=1,font="Helvetica",title="",outputfile="")
{
	links <- read.csv(filename)
	if(nrow(links)>0) { return }
	if ( chromosome > 0 ) links <- subset(links,chr==chromosome)

	#make sure these are numeric.
	links$start <- as.numeric(links$start)
	links$end <- as.numeric(links$end)
	links$cM <- as.numeric(links$cM)

	#build a list of our vertices
  people <- unique(rbind(unique(data.frame(name=links$u1,has.data=links$u1data,is.focal=links$u1focal)),unique(data.frame(name=links$u2,has.data=links$u2data,is.focal=links$u2focal))))

	g <- graph.data.frame(links,directed=FALSE,vertices=people)

  E(g)$width <- log(links$cM)
	if (use.edge.color) { E(g)$color <- links$chr }
	if (use.edge.label) { E(g)$label <- as.character(links$chr) }
	if (use.community)
	{
		com <- vertex.community(g)
		V(g)$color <- com$membership
	}
	else
	{
		V(g)$color <- people$has.data + 1 #add one so we don't get clear colors.
	}

	#handle title
	top.margin <- if (title!="") 1 else 0

	#now do the plot
	par(mai=c(0,0,top.margin,0))
	plot(g,vertex.label.cex=cex,edge.label.cex=cex,vertex.label.family=font,edge.label.family=font)

  title(title)
	if ( outputfile != "" ) dev.copy2pdf(file=outputfile)

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
