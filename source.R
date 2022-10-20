run_diffMap <- function(data=data, condition=condition, sigma="local"){
    destinyObj <- as.ExpressionSet(as.data.frame(t(data)))
    destinyObj$condition <- factor(condition)
    dm <- DiffusionMap(destinyObj, sigma)
    return(dm)
}

### 
plot_dm_3D <- function(dm=dm, dc=c(1:3), condition=condition, colour=colour){
    cond <- factor(condition)
    col <- factor(condition)
    levels(col) <- colour
    col <- as.vector(col)
    DCs <- paste("DC",dc, sep="")
    data <- data.frame(
        dm@eigenvectors[,DCs[1]], 
        dm@eigenvectors[,DCs[2]], 
        dm@eigenvectors[,DCs[3]]
    )
    colnames(data) <- DCs
    plot3d(data,col=col,alpha=0.8, size=5,box = FALSE)
	legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")
}

### 
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
		tmp <- rev(colorRampPalette(brewer.pal(11,'Spectral')[-6])( length(cells) ))
		names(tmp) <- cells[order(na.omit(slingPseudotime(crv)[,lineage]))]
		col <- tmp[ rownames(data) ]
		plot3d(data,col=col,alpha=0.8, size=5,box = FALSE)
	}else{
		condition <- condition[cells]
		col <- factor(condition)
		levels(col) <- colour
		col <- as.vector(col)
		plot3d(data,col=col,alpha=0.8, size=5,box = FALSE)
		legend3d("topright", legend = unique(condition), pch = 21, pt.bg = unique(col), cex=1.5, inset=c(0.02), bty = "n")
	}
}

### 
rankKeepNA <- function(x) {
    return(
        ifelse(
            is.na(x),
            NA,
            rank(
                x, 
                na.last = TRUE, 
                ties.method="random"
            )
        )
    )
}

get_pseudotime <- function(pseudotime, wthres=wthres){
    pseudoT <- list()
    for(lineage in 1:length(pseudotime@metadata$curves))local({
        curve <- pseudotime@metadata$curves[[lineage]]
        lambda <- curve$lambda
        weight <- curve$w
        ps <- curve$lambda
        ps[weight < wthres] <- NA
        ps <- rankKeepNA(ps)
        pseudoT[[lineage]] <<- ps
    })
    df <- t(do.call("rbind",pseudoT))
    colnames(df) <- names(pseudotime@metadata$curves)
    return(df)
}

###
plot_smoothed_gene_per_lineage <- function(
    rpkm_matrix=rpkm_matrix, 
    pseudotime=pseudotime, 
    lin=lin,
    geneName=geneName, 
    stages=stages, 
    clusters=clusters, 
    stage_colors=stage_colors,
    cluster_colors=cluster_colors,
    lineage_colors=lineage_colors
    ){

    pseudotime <- pseudotime[,lin]

    lineage_colors <- lineage_colors[lin]

    myplots <- list()
    total_pseudotime <- vector()

    if (length(lin)==1){

        lineage <- pseudotime
        total_pseudotime <- c(total_pseudotime, lineage)
        total_pseudotime <- na.omit(total_pseudotime)
        max_exp <- max(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,])
        max_pseudotime <- max(total_pseudotime)

        pseudotime_data <- data.frame(
            pseudotime=numeric(),
            lineage=numeric(),
            stages=character(),
            clusters=character(),
            gene=numeric()
        )


        lineage <- pseudotime
        sub_data <- rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]

        data <- data.frame(
            pseudotime=lineage,
            lineage=paste("Lineage ",lin, sep=""),
            stages=stages,
            clusters=clusters,
            gene=as.numeric(sub_data)
        )

        colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene")

        # print(lin)

        pseudotime_data <- rbind(pseudotime_data, data)

        p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
        geom_point(shape=21, size = 1.5,  aes(fill=clusters), color="white", na.rm = TRUE)+
        geom_point(shape=108, size = 6,  aes(y=-1, color=stages), na.rm = TRUE)+
        geom_smooth(color="black", na.rm = TRUE, method="loess", span=0.5)+
        ylab("log(CPM+1)") +
        theme_bw() +
        ggtitle(geneName)+
        # scale_color_manual(
        #     values=cluster_colors
        # ) +
        scale_fill_manual(
            values=cluster_colors
        ) +
        scale_color_manual(
            values=stage_colors
        ) +
        expand_limits(y = c(0,2))+
        theme(
            axis.text=element_text(size=16),
            axis.title=element_text(size=16),
            legend.text = element_text(size =16),
            legend.title=element_blank(),
            plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
            aspect.ratio=0.5,
            legend.position="bottom",
            strip.text.x = element_text(size = 16)
        )


    } else {

        for (lineages in 1:ncol(pseudotime)){
            lineage <- as.vector(pseudotime[,lineages])
            total_pseudotime <- c(total_pseudotime, lineage)
            total_pseudotime <- na.omit(total_pseudotime)
        }
    
        max_exp <- max(rpkm_matrix[rownames(rpkm_matrix) %in% geneName,])
        max_pseudotime <- max(total_pseudotime)

        pseudotime_data <- data.frame(
            pseudotime=numeric(),
            lineage=numeric(),
            stages=character(),
            clusters=character(),
            gene=numeric()
        )


        for (lineages in 1:ncol(as.data.frame(pseudotime))){
            lineage <- pseudotime[,lineages]
            sub_data <- rpkm_matrix[rownames(rpkm_matrix) %in% geneName,]

            data <- data.frame(
                pseudotime=lineage,
                lineage=paste("Lineage ",lineages, sep=""),
                stages=stages,
                clusters=clusters,
                gene=as.numeric(sub_data)
            )

            colnames(data) <- c("pseudotime", "lineage", "stages", "clusters", "gene")

            pseudotime_data <- rbind(pseudotime_data, data)
        }

        p <- ggplot(pseudotime_data, aes(x=pseudotime, y=gene))+
        geom_smooth(aes(group=lineage, color=lineage, fill=lineage), na.rm = TRUE, method="loess", span=0.5)+
        ylab("log(CPM+1)") +
        theme_bw() +
        ggtitle(geneName)+
        scale_color_manual(
            values=lineage_colors
        ) +
        scale_fill_manual(
            values=lineage_colors
        ) +
        expand_limits(y = c(0,2))+
        theme(
            axis.text=element_text(size=16),
            axis.title=element_text(size=16),
            legend.text = element_text(size =16),
            legend.title=element_blank(),
            plot.title = element_text(size=18, face="bold.italic", hjust = 0.5),
            aspect.ratio=0.5,
            legend.position="bottom",
            strip.text.x = element_text(size = 16)
        )
    }
    print(p)
}

smooth_gene_exp <- function(data=data, pseudotime=pseudotime, span=0.75){
    smooth_data <- data
    for (gene in 1:nrow(data)){
        gene_exp <- as.numeric(data[gene,])
        smooth <- loess(formula=gene_exp~pseudotime, span=span)
        smooth_data[gene,] <- predict(smooth, newdata=pseudotime)
    }
    return(smooth_data)
}
