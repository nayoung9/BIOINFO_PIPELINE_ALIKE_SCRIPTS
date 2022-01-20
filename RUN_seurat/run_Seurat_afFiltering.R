library(config)
library(Seurat)
library(dplyr)

arg <- commandArgs()
f_saved.rds <- arg[6]
rawData.SO <- readRDS(file=f_saved.rds)
new_param <- arg[7]
f_param <- new_param
config <- config::get(file=f_param)
config$dir.out <- arg[8]

if (length(config$dir.out) == 0){
	config$dir.out <- "./"
}else{
	cmd = paste("mkdir -p ", config$dir.out)
	system(cmd)
}

outFile <- function(x=""){
	str <- paste(config$dir.out,"/",x,sep="")
	return(str)
}
runStartTime <- gsub(" ","_",Sys.time())
print (paste("run start ",runStartTime, sep=""))
for (i in 1:length(config)){
	print(config[i])
}
convertpdf2png <- function(inpdf, outpng){
	cmd <- paste("convert -density 150 -trim -quality 100 -flatten -sharpen 0x1.0 -strip",inpdf, outpng)
    system(cmd)
}

#Filetering with UMI count / MT percent
#UMI upper/lower limits and MT percent limits
filtData.SO<- subset(rawData.SO, subset = nFeature_RNA > config$features.min & nFeature_RNA < config$features.max & percent.mt < config$percent.mt.max)
write.table(filtData.SO@meta.data, file=outFile("Matrix.afFilt.tsv"), quote=F, sep='\t', col.names=NA)
pdf(outFile("vlnplot_afFilt.pdf"), width=11, height=8)
VlnPlot(filtData.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
convertpdf2png(outFile("vlnplot_afFilt.pdf"),outFile("vlnplot_afFilt.png"))

print ("Normalization")
#normalization methods, and scale factor 
filtData.SO <- NormalizeData(filtData.SO, scale.factor = config$scale.factor)

#selection method, and variable counts 
print ("Find variable features")
Sys.setenv(R_CONFIG_ACTIVE="VariableGenes")
config <- config::get(file=f_param)
config$dir.out <- arg[8]
filtData.SO <- FindVariableFeatures(filtData.SO, mean.cutoff = c(config$mean.min,config$mean.max),dispersion.cutoff = c(config$dispersion.min, config$dispersion.max), num.bin= config$num.bin)
#print variable features
VariableFeatures <- filtData.SO@assays$RNA@meta.features
write.table(VariableFeatures, file=outFile("VariableFeatures.tsv"), sep='\t', col.names=NA, quote=F)
SortVariableFeatures <- subset(VariableFeatures, select = c("vst.variance.standardized", "vst.mean", "vst.variable"), subset=(vst.variable == TRUE))
write.table(SortVariableFeatures, file=outFile("OnlyVariableFeatures.tsv"), sep='\t', col.names=NA, quote=F)


#visualize features
pdf(outFile("VafiableFeatures.pdf"), width=11, height=8)
plot1 <- VariableFeaturePlot(filtData.SO)
#visualize strongly: variable features -> count n for parameter 
top10 <- head(VariableFeatures(filtData.SO),10)
plot2 <- LabelPoints(plot=plot1, points=top10, repel = TRUE)
print(plot1)
print(plot2)
dev.off()
convertpdf2png(outFile("VafiableFeatures.pdf"),outFile("VafiableFeatures.png"))

#scalining the data
all.genes <- rownames(filtData.SO)
filtData.SO <- ScaleData(filtData.SO, features = all.genes)


#PCA
Sys.setenv(R_CONFIG_ACTIVE="PCA")
config <- config::get(file=f_param)
config$dir.out <- arg[8]
print ("PCA")
filtData.SO <- RunPCA(filtData.SO, features=VariableFeatures(object=filtData.SO), npcs = config$npcs)
sink(outFile("PCA_result_top10.txt"))
print(filtData.SO[["pca"]], dims = 1:10, nfeatures = 10)
sink()
#visualize PCA Results 
pdf(outFile("PCA_results.pdf"), width=11, height=8);
VizDimLoadings(filtData.SO, dims=1, reduction="pca")
VizDimLoadings(filtData.SO, dims=2, reduction="pca")
DimPlot(filtData.SO, reduction = "pca")
for( i in 1:30){
	DimHeatmap(filtData.SO, dims=i, cells=500, balanced=TRUE)
}
DimHeatmap(filtData.SO, dims=1:30, cells=500, balanced=TRUE)
dev.off()
convertpdf2png(outFile("PCA_results.pdf"),outFile("PCA_results.png"))


#JaskStraw test: determine how many PCs to use 
filtData.SO <- JackStraw(filtData.SO, num.replicate = 100)
filtData.SO <- ScoreJackStraw(filtData.SO, dims=1:20)
pdf(outFile("JackStrawPlot.pdf"), width =11, height=8)
JackStrawPlot(filtData.SO, dims = 1:20)
dev.off()
convertpdf2png(outFile("JackStrawPlot.pdf"),outFile("JackStrawPlot.png"))

#draw elbowPlot
pdf(outFile("ElbowPlot.pdf"), width=11, height=8)
ElbowPlot(filtData.SO)
dev.off()
convertpdf2png(outFile("ElbowPlot.pdf"),outFile("ElbowPlot.png"))

##################################################################
# After Draw JacksStraw test, elbow plot 
# User should determine which PCs to use
# Are there any automation methods? 
# need to think about more.... 


print ("Clustering")
Sys.setenv(R_CONFIG_ACTIVE="Neighbors")
config <- config::get(file=f_param)
filtData.SO<- FindNeighbors(filtData.SO, dims =1:config$npcs ,k.param = config$K )
Sys.setenv(R_CONFIG_ACTIVE="Clustering")
config <- config::get(file=f_param)
config$dir.out <- arg[8]
# dims for clustering should be decided based on data: from before steps' results

#leiden
filtData.SO.leiden <- FindClusters(filtData.SO, resolution = config$resolution, algorithm=4)
write.table(filtData.SO.leiden@meta.data, outFile("/metadata.leiden.tsv"))

#louvain
filtData.SO.louvain <- FindClusters(filtData.SO, resolution = config$resolution, algorithm=1)
write.table(filtData.SO.louvain@meta.data, outFile("/metadata.louvain.tsv"))

Sys.setenv(R_CONFIG_ACTIVE="Embedding")
config <- config::get(file=f_param)
config$dir.out <- arg[8]
# non-linear dimension reduction ( for visualization only)
filtData.SO.leiden <- RunUMAP(filtData.SO.leiden, dims = 1:config$npcs)
pdf(outFile("UMAP_leiden.pdf"), width=11, height=8)
DimPlot(filtData.SO.leiden, reduction = "umap")
dev.off()
convertpdf2png(outFile("UMAP_leiden.pdf"),outFile("UMAP_leiden.png"))
filtData.SO.leiden <- RunTSNE(filtData.SO.leiden, dims = 1:config$npcs, perplexity = config$perplexity)
pdf(outFile("TSNE_leiden.pdf"), width=11, height=8)
DimPlot(filtData.SO.leiden, reduction = "tsne")
dev.off()
convertpdf2png(outFile("TSNE_leiden.pdf"),outFile("TSNE_leiden.png"))

# non-linear dimension reduction ( for visualization only)
filtData.SO.louvain <- RunUMAP(filtData.SO.louvain, dims = 1:config$npcs)
pdf(outFile("UMAP_louvain.pdf"), width=11, height=8)
DimPlot(filtData.SO.louvain, reduction = "umap")
dev.off()
convertpdf2png(outFile("UMAP_louvain.pdf"),outFile("UMAP_louvain.png"))
filtData.SO.louvain <- RunTSNE(filtData.SO.louvain, dims = 1:config$npcs,perplexity = config$perplexity)
pdf(outFile("TSNE_louvain.pdf"), width=11, height=8)
DimPlot(filtData.SO.louvain, reduction = "tsne")
dev.off()
convertpdf2png(outFile("TSNE_louvain.pdf"),outFile("TSNE_louvain.png"))

#Find markers per Clusters
Sys.setenv(R_CONFIG_ACTIVE="Marker_gene")
config <- config::get(file=f_param)
config$dir.out <- arg[8]
#nouse temp
if(FALSE){
filtData.SO.louvain <- FindAllMarkers(filtData.SO.louvain, only.pos = TRUE, min.pct = config$findMarker.min.pct, logfc.threshold = config$findMarker.logfc.thres)
write.table(filtData.SO.louvain %>% group_by(cluster), file=outFile("MarkerGenesPerClusters_louvain.txt"), quote=F, sep='\t', col.names=NA)
write.table(filtData.SO.louvain %>% group_by(cluster)%>% top_n(n = 5, wt = avg_logFC), file=outFile("MarkerGenesPerClusters_top5_louvain.txt"), quote=F, sep='\t', col.names=NA)
}

filtData.SO.leiden.MK <- FindAllMarkers(filtData.SO.leiden, only.pos = TRUE, logfc.threshold = config$fc.min, test.use=config$method)
write.table(filtData.SO.leiden.MK %>% group_by(cluster), file=outFile("MarkerGenesPerClusters_leiden.txt"), quote=F, sep='\t', col.names=NA)
write.table(filtData.SO.leiden.MK %>% group_by(cluster)%>% top_n(n = 5, wt = avg_logFC), file=outFile("MarkerGenesPerClusters_top5_leiden.txt"), quote=F, sep='\t', col.names=NA)

#top  
top5_percluster <- filtData.SO.leiden.MK %>% group_by(cluster)%>% top_n(n = 5, wt = avg_logFC)
top3_percluster <- filtData.SO.leiden.MK %>% group_by(cluster)%>% top_n(n = 3, wt = avg_logFC)


if(FALSE){
#Visualization
# This step also can be manually modifie
#
print ("Visualization")
pdf(outFile("vlnplots_perMarkerGenes_T3_all.pdf"), width=21, height=29.7)
VlnPlot(filtData.SO, features = c(top3_percluster$gene))
dev.off()
convertpdf2png(outFile("vlnplots_perMarkerGenes_T3_all.pdf"),outFile("vlnplots_perMarkerGenes_T3_all.png"))
pdf(outFile("vlnplots_perMarkerGenes_T3_count_all.pdf"), width=21, height=29.7)
VlnPlot(filtData.SO, features = c(top3_percluster$gene),slot = "counts", log = TRUE)
dev.off()
covertpdf2png(outFile("vlnplots_perMarkerGenes_T3_count_all.pdf"),outFile("vlnplots_perMarkerGenes_T3_count_all.png"))
pdf(outFile("FeaturePlot_perMarkerGenes_T3_all.pdf"), width=21, height=29.7)
FeaturePlot(filtData.SO, features = c(top3_percluster$gene))
dev.off()
convertpdf2png(outFile("FeaturePlot_perMarkerGenes_T3_all.pdf"),outFile("FeaturePlot_perMarkerGenes_T3_all.png"))

pdf(outFile("vlnplots_perMarkerGenes_T5_all.pdf"), width=21, height=29.7)
VlnPlot(filtData.SO, features = c(top5_percluster$gene))
dev.off()
convertpdf2png(outFile("vlnplots_perMarkerGenes_T5_all.pdf"),outFile("vlnplots_perMarkerGenes_T5_all.png"))
pdf(outFile("vlnplots_perMarkerGenes_T5_count_all.pdf"), width=21, height=29.7)
VlnPlot(filtData.SO, features = c(top5_percluster$gene),slot = "counts", log = TRUE)
dev.off()
convertpdf2png(outFile("vlnplots_perMarkerGenes_T5_count_all.pdf"),outFile("vlnplots_perMarkerGenes_T5_count_all.png"))
pdf(outFile("FeaturePlot_perMarkerGenes_T5_all.pdf"), width=21, height=29.7)
FeaturePlot(filtData.SO, features = c(top5_percluster$gene))
dev.off()
convertpdf2png(outFile("FeaturePlot_perMarkerGenes_T5_all.pdf"),outFile("FeaturePlot_perMarkerGenes_T5_all.png"))
}
#saveRDS(filtData.SO, file="after_filtered.rds")
save.image(file=outFile(paste(runStartTime,".RData",sep="")))

print ("Done")
