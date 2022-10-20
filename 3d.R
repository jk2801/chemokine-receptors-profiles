plot_dm_3D_lineage <- function(dm=dm, crv=crv, dc=c(1:3), lineage= "Lineage1", condition=condition, colour=colour, ifpesudo=FALSE){
    DCs <- paste("DC",dc, sep="")
    data <- data.frame(
        dm@eigenvectors[,DCs[1]], 
        dm@eigenvectors[,DCs[2]], 
        dm@eigenvectors[,DCs[3]]
    )
    colnames(data) <- DCs
	cells <- names(na.omit(slingPseudotime(crv)[,lineage]))
	data <- data[ cells, ]
	if( ifpesudo ){
		library(plotrix)
		tmp <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])( length(cells) ))
		names(tmp) <- cells[order(na.omit(slingPseudotime(crv)[,lineage]))]
		col <- tmp[ rownames(data) ]
		open3d(windowRect = c(10, 10, 500, 500))
		layout3d(matrix(1:2, 1,2), c(0.8, 0.2), 1)
		plot3d(data,col=col,alpha=0.8, size=5,box = FALSE)
		bgplot3d({
		plot.new()
		color.legend(1, 0.5, 1.05, 0.9, rect.col=tmp, legend=seq(0, 100, by=20), gradient="y", cex = 0.8)
  })
	}else{
		condition <- condition[cells]
		col <- factor(condition)
		levels(col) <- colour
		col <- as.vector(col)
		plot3d(data,col=col,alpha=0.8, size=5,box = FALSE)
		legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")
	}
}
