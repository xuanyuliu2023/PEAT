library(Seurat)
seuObj <- readRDS('human_all.rds')
seuObj_VAT <- subset(seuObj, depot == "VAT")
DefaultAssay(seuObj_VAT) <- 'RNA'

library(ggplot2)

# UMAP category data plot
categary='integrated_snn_res.0.4'
categary='integrated_snn_res.0.6'
categary="orig.ident"
categary="group"
categary="Phase"
categary="cell_type__ontology_label"
categary="cell_type__custom"
categary="cell_type2"
categary="cell_type"

categary.color.pallet  <- as.character(c(
  "#E82C61","#5B74F4",
  '#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c565B','#E236AF','#8C8C00','#700B35','#06BACE',
  "#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#4567EA"))
DimPlot(seuObj_VAT, label = T, repel = TRUE, reduction = "umap",group.by = categary, raster=F) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2')                        +
 scale_color_manual(values=categary.color.pallet)
ggsave(file=paste0(categary,'_UMAP_DimPlot.jpg'),dpi = 800,width=7.34,height=5.4,unit='in')
ggsave(file=paste0(categary,'_UMAP_DimPlot.pdf'))
save(seuObj_VAT,file='seuObj_VAT.Rdata')
# ------------------- reference mapping ------------
library(Seurat)

# change the current plan to access parallelization
plan("multiprocess", workers = 20)
#plan()
options(future.globals.maxSize = 20* 1000 * 1024^2)

library(ggplot2)

load('seuObj_VAT.Rdata')

seuObj.ref <- seuObj_VAT
seuObj.query <- get(load('/work/project/xuanyu/project/PEAT_V2/seurat/integrated/seuObj.integrated.ALL.RData'))
DefaultAssay(seuObj.ref) <- 'RNA'
seuObj.ref <- FindVariableFeatures(seuObj.ref,selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seuObj.query <- FindVariableFeatures(seuObj.query,selection.method = "vst", nfeatures = 3000, verbose = FALSE)
length(intersect(seuObj.ref@assays$RNA@var.features,seuObj.query@assays$RNA@var.features))

# Cell type classification using an integrated reference
anchors <- FindTransferAnchors(reference = seuObj.ref, query = seuObj.query,
    dims = 1:30, reference.reduction = "pca",normalization.method = 'LogNormalize')
predictions <- TransferData(anchorset = anchors, refdata = seuObj.ref$cell_type2,
    dims = 1:30)
seuObj.query <- AddMetaData(seuObj.query, metadata = predictions)

categary.color.pallet  <- as.character(c(
  '#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c564b','#E236AF','#8C8C00','#700B35','#06BACE',
"#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6367", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#bd0026","#1b9e77","#d95f02",
  "#7570b3","#e7298a","#e6ab02","#f5d206","#4a453f",
  "#69767c","#687f65","#e5e58b","#3caea3","#f6d55c"
))

celltype2_predict_info <- seuObj.query@meta.data$predicted.id
save(celltype2_predict_info,file='celltype2_predict_info.Rdata')

categary='predicted.id'
DimPlot(seuObj.query, label = T, repel = TRUE, reduction = "umap",group.by = categary) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2') +
 scale_color_manual(values=categary.color.pallet)
ggsave(file=paste0('celltype2','_ref_mapping.UMAP.jpg'),dpi = 800,width=5.6,height=4.8,unit='in')
ggsave(file=paste0('celltype2','_ref_mapping.UMAP_labelF.jpg'),dpi = 800,width=5.6,height=4.8,unit='in')
ggsave(file=paste0('celltype2','_ref_mapping.UMAP.pdf')

feature <- 'CYBB'
FeaturePlot(seuObj.query, features = feature, reduction = 'umap',cols =c('grey','red'))  +
 ggtitle(feature) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2')
ggsave(file=paste0(feature,'.UMAP.jpg'),width=5.6,height=4.82,unit='in',dpi = 600)


#projection of a query onto the reference UMAP structure

seuObj.query <- MapQuery(anchorset = seuObj.anchors, reference = seuObj.ref, query = seuObj.query,
    refdata = list(celltype = "subcluster"), reference.reduction = "pca", reduction.model = "umap")

save(seuObj.query,file='seuObj.query.Rdata')

# alluvial plot
df <- seuObj.query@meta.data[,c('cellType','predicted.id')]
library(dplyr)
df.wide <- df %>% group_by(cellType,predicted.id) %>% summarize(Freq=n())
library("ggalluvial")
ggplot(data = df.wide,
       aes(axis1 = cellType, axis2 = predicted.id,
           y = Freq)) +
  scale_x_discrete(limits = c("query", "reference"), expand = c(.2, .05)) +
  xlab("cellType") +
  geom_alluvium(aes(fill = cellType)) +
  geom_stratum(aes(fill = cellType)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() + scale_fill_manual(values=categary.color.pallet) +
  theme(legend.position="none")
ggsave(file="alluvial_plot.pdf")




# ------------------- reference mapping in a cell type ------------
library(Seurat)
load('../seuObj_VAT.Rdata')
# change the current plan to access parallelization
plan("multiprocess", workers = 20)
#plan()
options(future.globals.maxSize = 20* 1000 * 1024^2)

seuObj_subset <- subset(seuObj_VAT, cell_type2 == "t_cell")
seuObj.query <- get(load('/work/project/xuanyu/project/PEAT_V2/subclustering/TC/seuObj.integrated.TC.RData'))


library(ggplot2)
seuObj.ref <- seuObj_subset
DefaultAssay(seuObj.ref) <- 'RNA'
seuObj.ref <- FindVariableFeatures(seuObj.ref,selection.method = "vst", nfeatures = 3000, verbose = FALSE)

seuObj.query <- FindVariableFeatures(seuObj.query,selection.method = "vst", nfeatures = 3000, verbose = FALSE)
length(intersect(seuObj.ref@assays$RNA@var.features,seuObj.query@assays$RNA@var.features))

# Cell type classification using an integrated reference
anchors <- FindTransferAnchors(reference = seuObj.ref, query = seuObj.query,
    dims = 1:30, reference.reduction = "pca",normalization.method = 'LogNormalize')
predictions <- TransferData(anchorset = anchors, refdata = seuObj.ref$cell_type,
    dims = 1:30,k.weight=6)
seuObj.query <- AddMetaData(seuObj.query, metadata = predictions)

categary.color.pallet  <- as.character(c(
  '#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c564b','#E236AF','#8C8C00','#700B35','#06BACE',
"#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6367", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#bd0026","#1b9e77","#d95f02",
  "#7570b3","#e7298a","#e6ab02","#f5d206","#4a453f",
  "#69767c","#687f65","#e5e58b","#3caea3","#f6d55c"
))

celltype_predict_info <- seuObj.query@meta.data$predicted.id
save(celltype_predict_info,file='celltype_predict_info.Rdata')

categary='predicted.id'
DimPlot(seuObj.query, label = F, repel = TRUE, reduction = "umap",group.by = categary) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2') +
 scale_color_manual(values=categary.color.pallet)
ggsave(file='UMAP_DimPlot_combined.jpg',dpi=800)

DimPlot(object =seuObj.query, reduction = "umap", group.by = categary, split.by = categary, label = F) + 
theme_bw() + scale_color_manual(values=categary.color.pallet)
ggsave(file='UMAP_DimPlot_split_by_categary.jpg',dpi=800)

save(seuObj.query,file='seuObj.query.Rdata')

df <- seuObj.query@meta.data[,c('subcluster','predicted.id')]
library(dplyr)
df.wide <- df %>% group_by(subcluster,predicted.id) %>% summarize(Freq=n())
library("ggalluvial")
ggplot(data = df.wide,  aes(axis1 = subcluster, axis2 = predicted.id,  y = Freq)) +
  scale_x_discrete(limits = c("query", "reference"), expand = c(.2, .05)) +
  xlab("subcluster") +
  geom_alluvium(aes(fill = subcluster)) +
  geom_stratum(aes(fill = subcluster)) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() + scale_fill_manual(values=categary.color.pallet) +
  theme(legend.position="none")
ggsave(file="alluvial_plot.pdf")


