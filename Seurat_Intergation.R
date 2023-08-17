library(Seurat)
library(future)

# change the current plan to access parallelization
plan("multiprocess", workers = 10)
#plan()
options(future.globals.maxSize = 20* 1000 * 1024^2)

load('../QC/seuObj.list.30samps_scrublet.Rdata')
#sample names & aggr data path
sampleNames <- c("C11","C14","C15","C16","C01","C02","C04","C05","C06","C07","C08","S02","S03","S05","S08","S11","S12","S13","S14","S16","S17","S18","S19","S20","S21","S22","S23","S24","S25","S27")


# SCT transform
#seuObj.list <- lapply(X = seuObj.list, function(i){SCTransform(i,vars.to.regress = c("percent.mito"), 
method = "glmGamPoi",assay='RNA',vst.flavor = "v2")})
# Typically use 3,000 or more features for analysis downstream of sctransform
#features <- SelectIntegrationFeatures(object.list = seuObj.list, nfeatures = 3000)
# Run the PrepSCTIntegration() function prior to identifying anchors
#seuObj.list <- PrepSCTIntegration(object.list = seuObj.list, anchor.features = features)

#logNormalize the data
#summary(Matrix::colSums(seuObj.list[[1]]@assays$RNA@data))
for (i in 1:length(x = seuObj.list)) {
    seuObj.list[[i]] <- NormalizeData(object = seuObj.list[[i]],normalization.method = "LogNormalize",scale.factor = 10000)
}

#feature selection
for (i in 1:length(x = seuObj.list)) {
    seuObj.list[[i]] <- FindVariableFeatures(object = seuObj.list[[i]],
        selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    var.genes.ori <- seuObj.list[[i]]@assays$RNA@var.features
    print(length(x = VariableFeatures(object = seuObj.list[[i]])))
}

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = seuObj.list)
seuObj.list <- lapply(X = seuObj.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# integration
## ----------- CCA for small datasets
#cca.anchors <- FindIntegrationAnchors(object.list = seuObj.list, normalization.method = "SCT", anchor.features = features)
#seuObj <- IntegrateData(anchorset = cca.anchors, normalization.method = "SCT",new.assay.name = "integrated",dims = 1:30)

##############--------Fast integration using reciprocal PCA (RPCA) for large datasets
# find anchors;You can increase the strength of alignment by increasing the k.anchor parameter, which is set to 5 by default.
#Increasing this parameter to 20 will assist in aligning these populations.
rpca.anchors <- FindIntegrationAnchors(object.list = seuObj.list, anchor.features = features, reduction = "rpca", k.anchor = 5)
# integration
seuObj.integrated <- IntegrateData(anchorset = rpca.anchors,normalization.method = "LogNormalize",new.assay.name = "integrated",dims = 1:30)

#scale the data & Remove cell cycle effects
seuObj.integrated <- ScaleData(object = seuObj.integrated, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mito","S.Score", "G2M.Score"),
  display.progress = TRUE,model.use = "linear", do.scale = TRUE,do.center = TRUE)

#Perform dimensionality reduction by PCA and UMAP embedding
DefaultAssay(seuObj.integrated) <- "integrated"
seuObj.integrated <- RunPCA(seuObj.integrated, verbose = FALSE, assay = "integrated",npcs = 50)
seuObj.integrated <- RunUMAP(seuObj.integrated, reduction = "pca", assay = "integrated", dims = 1:30)
# clustering
seuObj.integrated <- FindNeighbors(seuObj.integrated, reduction = "pca", dims = 1:30)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 0.4)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 0.6)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 0.8)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 1.0)
# add meta info
library(plyr)

seuObj.integrated@meta.data$group <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sampleNames,
to=c("CTRL_old","CTRL_old","CTRL_old","CTRL_old","CTRL_old","CTRL_young","CTRL_young","CTRL_young","CTRL_old","CTRL_old","CTRL_old","CAD_without_T2D","CAD_without_T2D","CAD_without_T2D","CAD_with_T2D","CAD_with_T2D",
"CAD_without_T2D","CAD_without_T2D","CAD_with_T2D","CAD_with_T2D","CAD_without_T2D","CAD_without_T2D","CAD_with_T2D","CAD_with_T2D",
"CAD_with_T2D","CAD_without_T2D","CAD_with_T2D","CAD_with_T2D","CAD_with_T2D","CAD_with_T2D")
)


seuObj.integrated@meta.data$sex <- mapvalues(seuObj.integrated@meta.data$orig.ident,from=sampleNames,
to=c("male","male","male","male","female","male","male","male","male","male","female","male","female","male","male","male","male","male","male","male","male","male","male","male",
"male ","male","male","male","male","male")
)

# set Idents
Idents(seuObj.integrated) <- 'integrated_snn_res.0.6'

# Identify differential expressed genes across condition
# find all markers
# recorrect_umi=TRUE
markers.all <- FindAllMarkers(object =seuObj.integrated, assay = "SCT", only.pos = TRUE, min.pct = 0.25,
 logfc.threshold = 0.25, test.use="bimod", recorrect_umi=T)
write.table(markers.all,file="markers.all.recorrect_umi_T.tsv",sep="\t",row.names=F,quote=F,col.names=T)
# recorrect_umi=FALSE
markers.all <- FindAllMarkers(object =seuObj.integrated, assay = "SCT", only.pos = TRUE, min.pct = 0.25,
 logfc.threshold = 0.25, test.use="bimod", recorrect_umi=FALSE)
head(subset(markers.all,cluster == 'VSMC1'),n=20)
write.table(markers.all,file="markers.all.recorrect_umi_F.tsv",sep="\t",row.names=F,quote=F,col.names=T)

#Identify conserved cell type markers across conditions
Idents(seuObj.integrated) <- 'lineage'
adipocyte.markers <- FindConservedMarkers(seuObj.integrated, assay = "SCT", ident.1 = "Adipocyte", grouping.var = "group",verbose = FALSE)
head(adipocyte.markers)
# save
save(seuObj.integrated,file='seuObj.integrated.RData')

# visualization
library(ggplot2)

# UMAP category data plot
categary='integrated_snn_res.0.4'
categary='integrated_snn_res.0.6'
categary="orig.ident"
categary="group"
categary="Phase"
categary.color.pallet  <- as.character(c(
  "#E82C61","#5B74F4",
  '#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c565B','#E236AF','#8C8C00','#700B35','#06BACE',
  "#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#4567EA"))
DimPlot(seuObj.integrated, label = T, repel = TRUE, reduction = "umap",group.by = categary) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2')                        +
 scale_color_manual(values=categary.color.pallet)
ggsave(file=paste0(categary,'_UMAP_DimPlot.jpg'),width=5.6,height=4.8,unit='in',dpi = 800)

DimPlot(object =seuObj.integrated, reduction = "umap", group.by = "subgroup",split.by = "subgroup", label = F) + theme_bw() + scale_color_manual(values=cluster.color.pallet)
ggsave(file='UMAP_DimPlot_split_by_subgroup.jpg',width=14.6,height=3.92,units='in',dpi=800)
ggsave(file='UMAP_DimPlot_split_by_subgroup.pdf')


# feature UMAP plot
feature='LPL'
feature='CD14'

feature <- 'percent.mito'
feature <- 'nFeature_RNA'
feature <- 'nCount_RNA'
feature <- 'percent.ribo'

FeaturePlot(seuObj.integrated, features = feature, reduction = 'umap',cols =c('grey','red'),raster=T)  +
 ggtitle(feature) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2')
ggsave(file=paste0(feature,'.UMAP.jpg'),width=5.6,height=4.82,unit='in',dpi = 600)


# violin plot
#stack violin plot
features = c('LPL','ADIPOQ','PPARG','DCN','PDGFRA','MYH11','KCNJ8','VWF','EMCN','MSLN',
'KRT19','NRXN1','PTPRC','CD14','CD68','IL7R','FLT3','MS4A1','KLRD1','CPA3')
# best markers
features= c('LPL','ADIPOQ',
'DCN','PDGFRA',
'CD14','CD163',
'CPA3',
'FLT3',
'MCR1','F13A1',
'MS4A1','BANK1',
'SKAP1','IL7R','KLRD1',
'VWF','EMCN',
'PROX1','MMRN1',
'MSLN','KRT19',
'PTPRC'
)

violin.color.pallet  <- as.character(c(
  "#E82C61","#5B74F4",
  '#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c565B','#E236AF','#8C8C00','#700B35','#06BACE',
  "#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#4567EA"))

VlnPlot(seuObj.integrated, assay = 'SCT', features = features,stack= TRUE,flip=TRUE) + theme_bw() +
theme(legend.position='none') +  scale_fill_manual(values=violin.color.pallet) + labs(x='',y='Normalized gene expression')
ggsave(file=paste0('modulatedSMC','.feature.violinPlot.pdf'))
ggsave(file=paste0('marker','.feature.violinPlot.pdf'))

# cells to be discarded
cellsTobeDiscarded <- rownames(subset(seuObj.integrated@meta.data, integrated_snn_res.0.6 %in% c('8','14','25','22','29','28','23')))
save(cellsTobeDiscarded,file='cellsTobeDiscarded.RData')

