counts <- Read10X(data.dir =dir)
sce <- CreateSeuratObject(counts,min.cells = 100, min.features = 300)
table(sce@meta.data$orig.ident)  
sce[[]]
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HBm <- match(HB.genes, rownames(sce@assays$RNA)) 
HB.genes <- rownames(sce@assays$RNA)[HBm] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
sce[["percent.HB"]]<-PercentageFeatureSet(sce, features=HB.genes)
col.num <- length(levels(sce@active.ident))
violin <- VlnPlot(sce,
                  features = c("nFeatureRNA", "nCountRNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, #
                  ncol = 4) + 
  theme(axis.title.x=elementblank(), axis.text.x=elementblank(), axis.ticks.x=elementblank())

violin
ggsave("vlnplotbeforeqc.pdf", plot = violin, width = 12, height = 6) 
ggsave("vlnplotbeforeqc.png", plot = violin, width = 12, height = 6)  

plot1=FeatureScatter(sce, feature1 = "nCountRNA", feature2 = "percent.mt")
plot2=FeatureScatter(sce, feature1 = "nCountRNA", feature2 = "nFeatureRNA")
plot3=FeatureScatter(sce, feature1 = "nCountRNA", feature2 = "percent.HB")
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
plot1

plot2

plot3

pearplot

sce <- subset(sce, subset = nFeatureRNA > 300& nFeatureRNA < 7500 & percent.mt < 10 & percent.HB < 3 & nCountRNA < 100000)
sce
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)


sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000) 
# Identify the 10 most highly variable genes，
top10 <- head(VariableFeatures(sce), 10) 
# plot variable features with and without labels  
plot1 <- VariableFeaturePlot(sce) 

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 

plot 


scale.genes <-  rownames(sce)
sce <- ScaleData(sce, features = scale.genes)

scale.genes <-  VariableFeatures(sce)
sce <- ScaleData(sce, features = scale.genes)

#GetAssayData(sce,slot="counts",assay="RNA") 
              
#GetAssayData(sce,slot="data",assay="RNA")

#GetAssayData(sce,slot="scale.data",assay="RNA") 


sce <- RunPCA(sce, features = VariableFeatures(sce)) 
plot1 <- DimPlot(sce, reduction = "pca", group.by="orig.ident") 

plot1 
 Determine the ‘dimensionality’ of the dataset 
###ElbowPlot() 
plot2 <- ElbowPlot(sce, ndims=50, reduction="pca") 

plot2

plotc <- plot1+plot2
ggsave("pca.pdf", plot = plotc, width = 8, height = 4) 
ggsave("pca.png", plot = plotc, width = 8, height = 4)

pc.num=1:50


###Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm. First calculate k-nearest neighbors and construct the SNN graph. Then optimize the modularity function to determine clusters. For a full description of the algorithms, see Waltman and van Eck (2013) The European Physical Journal B. Thanks to Nigel Delaney (evolvedmicrobe@github) for the rewrite of the Java modularity optimizer code in Rcpp!
sce <- FindNeighbors(sce, dims = pc.num) 
###
sce <- FindClusters(sce, resolution = 0.8)

table(sce@meta.data$seurat_clusters)
metadata <- sce@meta.data
cellcluster <- data.frame(cellID=rownames(metadata), clusterID=metadata$seuratclusters)
sce = RunTSNE(sce, dims = pc.num)
embedtsne <- Embeddings(sce, 'tsne')
write.csv(embedtsne,'embedtsne.csv')
plot1 = DimPlot(sce, reduction = "tsne") 
##
plot1
#
DimPlot(sce, reduction = "tsne",label = TRUE) 

ggsave("tSNE.pdf", plot = plot1, width = 8, height = 7)


#
sce <- RunUMAP(sce, dims = pc.num)
embedumap <- Embeddings(sce, 'umap')
write.csv(embedumap,'embedumap.csv') 
plot2 = DimPlot(sce, reduction = "umap") 
plot2
ggsave("UMAP.pdf", plot = plot2, width = 8, height = 7)


plotc <- plot1+plot2+ plotlayout(guides = 'collect')
plotc
ggsave("tSNEUMAP.pdf", plot = plotc, width = 10, height = 5)

diff.wilcox = FindAllMarkers(sce)

all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(pval<0.05)
top10 = all.markers %>% groupby(cluster) %>% topn(n = 10, wt = avglogFC)

sce@meta.data$cellType[sce$seuratclusters%in%c('15')]<-"B cells"
