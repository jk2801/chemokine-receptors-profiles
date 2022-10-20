marker.gene.set=read.xlsx("gene_set.xlsx",sheetIndex = 1,header = F)
colnames(marker.gene.set)=c("gene","set")
marker.gene.set=marker.gene.set[marker.gene.set$gene %in% rownames(T_cell),]

for (i in unique(marker.gene.set$set)) {
  marker.gene.set_small=marker.gene.set%>%filter(set==i)
  genes.for.scoring <- list(marker.gene.set_small$gene)
  T_cell <- AddModuleScore(object = T_cell, features = genes.for.scoring, name = i)
}
ggplot(tsne_data, aes(x=tSNE_1,y=tSNE_2))+
  geom_density_2d_filled(bins=10)+
  scale_fill_manual(values = colorRampPalette(c("#ffffff","#75aadb"))(10))+
  geom_point(aes(color=ident),alpha=0.5,size=0.4,shape=16)+
  scale_color_manual(values = color_sample)+
  scale_x_continuous(limits=c(-33,32),expand = c(0.01,0))+
  scale_y_continuous(limits=c(-34.5,30.5),expand = c(0.01,0))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
  ggplot(umap_data, aes(x=UMAP_1,y=UMAP_2))+
  geom_density_2d_filled(bins=10)+
  scale_fill_manual(values = colorRampPalette(c("#ffffff","#75aadb"))(10))+
  geom_point(aes(color=ident),alpha=0.5,size=0.4,shape=16)+
  scale_color_manual(values = color_cluster)+
  scale_x_continuous(expand = c(0.01,0))+
  scale_y_continuous(expand = c(0.01,0))+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )
  FeaturePlot(T_cell,features = c("Naive_markers1","Inhibitory_signature1","Effector_molecular1","Co_stimulatory_moleculars1","Transcription_factors1","Treg_marekers1"),
            reduction = "tsne",pt.size = 0.5,cols = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"),order = T)&
  theme_bw()&
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5,size = 10)
  )
  library(RColorBrewer)
library(scales)

bubble_gene=read.table("./fig1/bubble_gene.txt",header = F,sep = "\t",stringsAsFactors = F)
colnames(bubble_gene)=c("cluster","gene")
bubble_gene$cluster=factor(bubble_gene$cluster,levels = c("B cell","CAFs","Epithelial cell","Myeloid","T cell"))
bubble_gene=bubble_gene%>%arrange(cluster,gene)

DotPlot(T_cell,features = bubble_gene$gene)+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))+
  coord_flip()+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
  VlnPlot(T_cell,features = bubble_gene$gene,pt.size = 0,stack = T)&
  scale_fill_manual(values = c(brewer.pal(11,"Set3"),brewer.pal(3,"Dark2")))&
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    strip.text.x = element_text(angle = 45,size = 10,hjust = 0,vjust = 0),
    legend.position = "none"
  )
