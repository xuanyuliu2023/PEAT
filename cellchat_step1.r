library(plyr)
library(CellChat)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
options(stringsAsFactors = FALSE)
library(Seurat)

# Prepare input data for CelChat analysis
load('/work/project/xuanyu/project/PEAT_V2/seurat/integrated/seuObj.integrated.ALL.RData')
> table(seuObj.integrated@meta.data$cellType)

 Macrophage   Adipocyte        ASPC          TC          BC         VEC
      17366       15893       15837        6682        3204        3062
        LEC      Neural Mesothelial    Pericyte          NK      DC_Neu
       2484        2198        2100        1660         844         726
     Plasma         SMC        Mast
        652         358         320

#meta.data <- subset(seuObj.integrated@meta.data,cellType %in% c('vEDO','FB','CM','Peri','macrophages','SMC','Tcell','Neuronal'))
meta.data <- seuObj.integrated@meta.data

> table(meta.data$cellType)[order(table(meta.data$cellType),decreasing=T)]

 Macrophage   Adipocyte        ASPC          TC          BC         VEC
      17366       15893       15837        6682        3204        3062
        LEC      Neural Mesothelial    Pericyte          NK      DC_Neu
       2484        2198        2100        1660         844         726
     Plasma         SMC        Mast
        652         358         320

> table(meta.data$cellType,meta.data$group)

              CTRL CAD_without_T2D CAD_with_T2D
  Macrophage  3694            7095         6577
  Adipocyte   4734            4795         6364
  ASPC        4686            5517         5634
  TC          3098            1692         1892
  BC          2295             308          601
  VEC         1382             817          863
  LEC         1432             529          523
  Neural       197             787         1214
  Mesothelial 1016             437          647
  Pericyte     590             474          596
  NK           313             255          276
  DC_Neu       134             283          309
  Plasma       372             117          163
  SMC          118              93          147
  Mast          83             139           98

allGenes <- rownames(seuObj.integrated)
geneTypeInfo <- read.table(file='/work/project/xuanyu/resource/10XGenomics/refdata-gex-GRCh38-2020-A/genes.gtf.geneType.tsv',header=F,sep='\t',stringsAsFactors=F)
protein_coding_genes <- subset(geneTypeInfo,V2=="protein_coding")$V1
allProteinCodingGenes <- allGenes[allGenes %in% protein_coding_genes]

meta.data$anno <- meta.data$cellType

meta.CAD_without_T2D <- subset(meta.data,group=="CAD_without_T2D")
meta.CAD_without_T2D <- meta.CAD_without_T2D[,c('orig.ident','cellType','group','anno')]
colnames(meta.CAD_without_T2D) <- c('orig.ident','labels','condition','anno')
cells.CAD_without_T2D <- rownames(meta.CAD_without_T2D)
data.input.CAD_without_T2D <- seuObj.integrated@assays$RNA@data[allProteinCodingGenes,cells.CAD_without_T2D]
unique(meta.CAD_without_T2D$labels) # check the cell labels
factorOrder <- c("Macrophage", "Adipocyte", "ASPC", "TC", "BC", "VEC", "LEC", "Neural", "Mesothelial", "Pericyte", "NK", "DC_Neu", "Plasma", "SMC", "Mast")
meta.CAD_without_T2D$anno <- factor(meta.CAD_without_T2D$anno,levels=factorOrder)
meta.CAD_without_T2D$labels <- factor(meta.CAD_without_T2D$labels,levels=factorOrder)
dim(meta.CAD_without_T2D)

meta.CAD_with_T2D <- subset(meta.data,group=="CAD_with_T2D")
meta.CAD_with_T2D <- meta.CAD_with_T2D[,c('orig.ident','cellType','group','anno')]
colnames(meta.CAD_with_T2D) <- c('orig.ident','labels','condition','anno')
cells.CAD_with_T2D <- rownames(meta.CAD_with_T2D)
data.input.CAD_with_T2D <- seuObj.integrated@assays$RNA@data[allProteinCodingGenes,cells.CAD_with_T2D]
unique(meta.CAD_with_T2D$labels) # check the cell labels
factorOrder <- c("Macrophage", "Adipocyte", "ASPC", "TC", "BC", "VEC", "LEC", "Neural", "Mesothelial", "Pericyte", "NK", "DC_Neu", "Plasma", "SMC", "Mast")
meta.CAD_with_T2D$anno <- factor(meta.CAD_with_T2D$anno,levels=factorOrder)
meta.CAD_with_T2D$labels <- factor(meta.CAD_with_T2D$labels,levels=factorOrder)
dim(meta.CAD_with_T2D)


meta.CTRL <- subset(meta.data,group=="CTRL")
meta.CTRL <- meta.CTRL[,c('orig.ident','cellType','group','anno')]
colnames(meta.CTRL) <- c('orig.ident','labels','condition','anno')
cells.CTRL <- rownames(meta.CTRL)
data.input.CTRL <- seuObj.integrated@assays$RNA@data[allProteinCodingGenes,cells.CTRL]
unique(meta.CTRL$labels) # check the cell labels
factorOrder <- c("Macrophage", "Adipocyte", "ASPC", "TC", "BC", "VEC", "LEC", "Neural", "Mesothelial", "Pericyte", "NK", "DC_Neu", "Plasma", "SMC", "Mast")
meta.CTRL$anno <- factor(meta.CTRL$anno,levels=factorOrder)
meta.CTRL$labels <- factor(meta.CTRL$labels,levels=factorOrder)
dim(meta.CTRL)

# Create a CellChat object
#If input is a Seurat or SingleCellExperiment object, the meta data in the object will be used by default and USER must provide group.by to define the cell groups. e.g, group.by = “ident” for the default cell identities in Seurat object.
cellchat.CAD_without_T2D <- createCellChat(object = data.input.CAD_without_T2D, meta = meta.CAD_without_T2D, group.by = "labels")
cellchat.CAD_with_T2D <- createCellChat(object = data.input.CAD_with_T2D, meta = meta.CAD_with_T2D, group.by = "labels")
cellchat.CTRL <- createCellChat(object = data.input.CTRL, meta = meta.CTRL, group.by = "labels")


save(cellchat.CAD_without_T2D,file="cellchat.CAD_without_T2D.ori.Rdata")
save(cellchat.CAD_with_T2D,file="cellchat.CAD_with_T2D.ori.Rdata")
save(cellchat.CTRL,file="cellchat.CTRL.ori.Rdata")
