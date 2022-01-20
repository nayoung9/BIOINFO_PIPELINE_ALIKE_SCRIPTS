library(config)
library(Seurat)
library(plyr)

arg <- commandArgs();
dir_input = arg[6] #basically 10X 
f_param = arg[7]
config <- config::get(file = f_param)

cdir <- getwd()
cmd = paste("mkdir -p  ",cdir,"/SEURAT/BF_filtering/",sep="")
system(cmd)

outFile<- function(x=""){
	str <- paste(cdir,"/SEURAT/BF_filtering/",x,sep="")
	return (str)
}
convertpdf2png <- function(inpdf, outpng){
	cmd <- paste("convert -density 150 -trim -quality 100 -flatten -sharpen 0x1.0 -strip ",inpdf," " ,outpng)
	system(cmd)
}

print(paste("READ: ",dir_input,sep=""))
#argument 1: input 
# 10X input: directory 
# ... searching for additional available input 
#argument 2: basename of output files(projectname)

rawData <- Read10X(data.dir = dir_input)
#min filtering parameters:min.cells, min.features
rawData.SO <- CreateSeuratObject(counts=rawData, min.cells=config$cells.min, min.features=config$features.min) 

rawData.SO<- subset(rawData.SO, subset = nFeature_RNA > config$features.min & nFeature_RNA < config$features.max)
#calc MT data
rawData.SO[["percent.mt"]] <- PercentageFeatureSet(rawData.SO, pattern ="^MT-")
write.table(rawData.SO@meta.data, file=outFile("Matrix.bfFilt.tsv"), quote=F, sep='\t', col.names=NA)
pdf(outFile("vlnplot_bfFilt.pdf"), width=11, height=8)
VlnPlot(rawData.SO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
convertpdf2png(outFile("vlnplot_bfFilt.pdf"),outFile("vlnplot_bfFilt.png"))

#print stat file 
sink(outFile("BF_stat.txt"))
paste("Total feature count: ", lengths(rawData.SO@assays$RNA@data@Dimnames[1]), sep="")
paste("Total cell count: ", lengths(rawData.SO@assays$RNA@data@Dimnames[2]), sep="")
cat("Summary of MT ratio")
summary(rawData.SO[["percent.mt"]])
cat("Summary of feature counts")
summary(rawData.SO[["nFeature_RNA"]])
sink()


saveRDS(rawData.SO, file=outFile("raw_data.rds"))
