load('/work/project/xuanyu/project/PEAT_latest/seurat/final_round3/meta.data.table.RData')
sampleNames <- c("C11","C14","C15","C16","C01","C02","C04","C05","C06","C07","C08","S02","S03","S05",
"S08","S11","S12","S13","S14","S16","S17","S18","S19","S20","S21","S22","S23","S24","S25","S27")

load('/work/project/xuanyu/project/PEAT_latest/seurat/QC/seuObj.list.30samps_scrublet.Rdata')
meta.df<- subset(meta.data.table,cellType =="Adipocyte")[,c('cellID','subcluster')]
selectedCells<- rownames(subset(meta.data.table,cellType =="Adipocyte"))
load('cell2discard_59.Rdata')
cellsToKeep <- selectedCells[!selectedCells %in% cell2discard_59]


library(Seurat)
# change the current plan to access parallelization
plan("multiprocess", workers = 10)
#plan()
options(future.globals.maxSize = 20* 1000 * 1024^2)

# filter cells
seuObj.list <- lapply(seuObj.list,function(x){AddMetaData(x, row.names(x@meta.data), col.name = 'cellID')})
seuObj.list <- lapply(seuObj.list,function(x){subset(x, cellID %in% cellsToKeep)})


## ----------- harmony integration for small datasets
library(harmony)
library(scCustomize)
seuObj.integrated <- Merge_Seurat_List(
  seuObj.list,
  add.cell.ids = NULL,
  merge.data = TRUE,
  project = "SeuratProject"
)

seuObj.integrated <- NormalizeData(seuObj.integrated,normalization.method = "LogNormalize",
    scale.factor = 10000)
seuObj.integrated <- FindVariableFeatures(seuObj.integrated,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seuObj.integrated <- ScaleData(object = seuObj.integrated, vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mito","S.Score", "G2M.Score"),
  display.progress = TRUE,model.use = "linear", do.scale = TRUE,do.center = TRUE)

seuObj.integrated <- RunPCA(seuObj.integrated, verbose = FALSE,npcs = 50)
seuObj.integrated <- RunHarmony(seuObj.integrated, group.by.vars="orig.ident")
seuObj.integrated <- RunUMAP(seuObj.integrated , reduction = "harmony", dims = 1:20)

# clustering
seuObj.integrated <- FindNeighbors(seuObj.integrated, reduction = "harmony", dims = 1:20)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 0.2)
seuObj.integrated <- FindClusters(seuObj.integrated, resolution = 0.3)
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
"male","male","male","male","male","male")
)

seuObj.integrated@meta.data$group <- factor(seuObj.integrated@meta.data$group,levels=c('CTRL_young','CTRL_old','CAD_without_T2D','CAD_with_T2D'))

seuObj.integrated@meta.data$subcluster <- paste0('Ad',seuObj.integrated@meta.data$RNA_snn_res.0.4)
seuObj.integrated@meta.data$subcluster <- factor(seuObj.integrated@meta.data$subcluster,levels=c('Ad0','Ad1','Ad2','Ad3','Ad4'))
# set Idents
Idents(seuObj.integrated) <- 'subcluster'


# Phylogenetic Analysis of Identity Classes
library(ape)
seuObj.integrated <- BuildClusterTree(object = seuObj.integrated,assay='RNA',
features=row.names(seuObj.integrated@assays$RNA))
plot(Tool(object = seuObj.integrated, slot = 'BuildClusterTree'))

# Identify differential expressed genes across condition
markers.all <- FindAllMarkers(object =seuObj.integrated, assay = "RNA", only.pos = TRUE, min.pct = 0.25,
 logfc.threshold = 0.25, test.use="bimod")
write.table(markers.all,file="markers.all.tsv",sep="\t",row.names=F,quote=F,col.names=T)
head(subset(markers.all,cluster == 'Ad3'),n=20)

#Identify conserved cell type markers across conditions
Idents(seuObj.integrated) <- 'lineage'
adipocyte.markers <- FindConservedMarkers(seuObj.integrated, assay = "SCT", ident.1 = "Adipocyte", grouping.var = "group",verbose = FALSE)
head(adipocyte.markers)


# save
save(seuObj.integrated,file='seuObj.integrated.Adipocyte.RData')

# visualization
library(ggplot2)

# UMAP category data plot
categary='RNA_snn_res.0.2'
categary='RNA_snn_res.0.4'
categary="orig.ident"
categary="group"
categary="Phase"
categary="subcluster"

categary.color.pallet  <- as.character(c(
  '#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c565B','#E236AF','#8C8C00','#700B35','#06BACE',
  "#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#4567EA", "#E82C61","#5B74F4"))

DimPlot(seuObj.integrated, label = F, raster=F, repel = TRUE, reduction = "umap",group.by = categary) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x='UMAP 1', y='UMAP 2')                        +
 scale_color_manual(values=categary.color.pallet)
ggsave(file=paste0(categary,'_UMAP_DimPlot.jpg'),width=4.33,height=3.7,unit='in',dpi = 800)

cluster.color.pallet <- c('#1f77b4','#ff7f0e','#2ca02c','#d62728')
DimPlot(object =seuObj.integrated, reduction = "umap", group.by = "group",split.by = "subgroup", label = F) + theme_bw() + scale_color_manual(values=cluster.color.pallet)
ggsave(file='UMAP_DimPlot_split_by_group.jpg',width=14.6,height=3.92,units='in',dpi=800)



# feature UMAP plot
feature='UCP1'
feature='DDX58'

feature <- 'percent.mito'
feature <- 'nFeature_RNA'
feature <- 'nCount_RNA'
feature <- 'percent.ribo'

FeaturePlot(seuObj.integrated, features = feature, reduction = 'umap',cols =c('grey','red'),raster=T)  +
 ggtitle(feature) + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x='UMAP 1', y='UMAP 2')
ggsave(file=paste0(feature,'.UMAP.jpg'),width=5.6,height=4.82,unit='in',dpi = 600)

# density_plot_for_sparse_marker
library(Nebulosa)
feature <- 'percent.mito'
feature <- 'nFeature_RNA'
feature <- 'nCount_RNA'
feature <- 'percent.ribo'

feature='DGAT2'
plot_density(seuObj.integrated, feature, pal="magma",size = 1) +theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   labs(x='UMAP 1', y='UMAP 2')
ggsave(file=paste0(feature,'.marker_density_plot.pdf'),width=5.22,height=3.89)
ggsave(file=paste0(feature,'.marker_density_plot.jpg'),width=5.22,height=3.89,unit='in',dpi = 600)

featureVector <- c('DGAT2','PDE5A','FBLN1','FOSB','PTCHD4')
featureVector <- c('percent.mito','nFeature_RNA','nCount_RNA','percent.ribo')
for (theFeature in featureVector) {
print(theFeature)
plot_density(seuObj.integrated, theFeature, pal="magma",size = 0.6) +theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   labs(x='UMAP 1', y='UMAP 2')
ggsave(file=paste0(theFeature,'.marker_density_plot.jpg'),width=5.6,height=4.82,unit='in',dpi = 600)
ggsave(file=paste0(theFeature,'.marker_density_plot.pdf'),width=4.54,height=3.72)      
}


# violin plot
#stack violin plot
features = c('LPL','ADIPOQ','PPARG','DCN','PDGFRA','MYH11','KCNJ8','VWF','EMCN','MSLN',
'KRT19','NRXN1','PTPRC','CD14','CD68','IL7R','FLT3','MS4A1','KLRD1','CPA3d
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
'PTPRC',
'STEAP4', 'MYH11',
'NRXN1'
)

features= c('UCP1','PRDM16','ADIPOQ')
features=c('ADIPOQ','GALNT13','TNFSF10','PNPLA3','GRIA4','PGAP1','EBF2','AGMO')
features=c('ADIPOQ','PTPRC','VWF','PDGFRA')
features=c('LPL','PDE3B','AQP7','ABCA9','ABCA6','NOD1','MYCBP2','DDB2','ZMAT3')
features=c('UCP1','PRDM16','PPARGC1A','CIDEA','SLC6A8','GATM','GAMT','CKMT1A','CKMT1B','CKMT2','CKB','EBR2','PPARGC1A','ESRRG')
features=c('ADIPOQ','PDGFRA','VWF','PROX1','MYH11','SDC1','CD19','MSLN','STEAP4','RBFOX3','CD3E','KLRD1','LYVE1','CPA3','CD1C')
violin.color.pallet  <- as.character(c(
  "#E82C61","#5B74F4",
  '#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c565B','#E236AF','#8C8C00','#700B35','#06BACE',
  "#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#4567EA"))

VlnPlot(seuObj.integrated, assay = 'RNA', features = features,stack= TRUE,flip=TRUE) + theme_bw() +
theme(legend.position='none') +  scale_fill_manual(values=violin.color.pallet) +
labs(x='',y='Normalized gene expression')

ggsave(file=paste0('marker','.feature.violinPlot.pdf'))

# DotPlot
Features <- c("PRSS23", "PDE5A", "CAV2", "GLRX", "NTRK3", "DGAT2", "MME", "CLMN", "MYADM", "RFX2", "FOSB", "LMOD1", "FBLN1", "RCAN1", "TAGLN", "TPM4", "KCNMB2-AS1", "PTCHD4", "MIR34AHG", "ZMAT3")
DotPlot(seuObj.integrated, assay = 'RNA', features = Features,cols = c("lightgrey", "red"))  + theme_bw() + 
labs(x='',y='Normalized gene expression') + 
theme(axis.text.x = element_text(angle = 30, size = 6, hjust = 0.5, vjust = 0.5))
ggsave(file=paste0('subcluster_marker','.feature.DotPlot.pdf'))


# ----------------meta.data--------------
# prepare meta data table
emb.umap <- FetchData(object = seuObj.integrated, vars = c("UMAP_1", "UMAP_2"))
meta.data <- seuObj.integrated@meta.data
meta.data$cellID <- row.names(meta.data)
meta.data.table <- cbind(emb.umap,meta.data)
save(meta.data.table,file="meta.data.table.RData")

NormData <- seuObj.integrated@assays$RNA@data
# split the sparse matrix to avoid the error
RNAdata <- NormData
halfColNum <- floor(ncol(RNAdata)/2)
RNAdata1 <- RNAdata[,1:halfColNum]
RNAdata2 <- RNAdata[,(halfColNum+1):ncol(RNAdata)]
raw.X.data1 <- t(as.matrix(RNAdata1))
raw.X.data2 <- t(as.matrix(RNAdata2))
raw.X.data <- rbind(raw.X.data1,raw.X.data2)
meta.data.NormDat <- cbind(meta.data.table,raw.X.data)
save(meta.data.NormDat, file='meta.data.NormDat.Rdata')

DefaultAssay(seuObj.integrated) <- 'RNA'
#split Dim group plot
DimPlot(object = seuObj.integrated, raster=FALSE, reduction = "umap", group.by = "sex",label = F,split.by ="sex") + theme_bw() + 
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values=c( '#E82C61','#5B74F4'))
ggsave(file='split_DimPlot_bySex.jpg')
ggsave(file='split_DimPlot_bySex.pdf')
DimPlot(object = seuObj.integrated, reduction = "umap", raster=FALSE, group.by = "group",label = F,split.by ="group") + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values=c('#1f77b4','#ff7f0e','#2ca02c','#d62728'))
ggsave(file='split_DimPlot_by4groups.jpg')
ggsave(file='split_DimPlot_by4groups.pdf')
DimPlot(object = seuObj.integrated, reduction = "umap", raster=FALSE, group.by = "Phase",label = F,split.by ="Phase") + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values=c('#1f77b4','#ff7f0e','#2ca02c'))
ggsave(file='split_DimPlot_byPhase.jpg')
ggsave(file='split_DimPlot_byPhase.pdf')

#join Dim group plot
DimPlot(object = seuObj.integrated, raster=FALSE, reduction = "umap", group.by = "sex",label = F) + theme_bw() + 
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values=c( '#E82C61','#5B74F4'))
ggsave(file='join_DimPlot_bySex.jpg')
ggsave(file='join_DimPlot_bySex.pdf')
DimPlot(object = seuObj.integrated, reduction = "umap", raster=FALSE, group.by = "group",label = F) + theme_bw() +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values=c('#1f77b4','#ff7f0e','#2ca02c','#d62728'))
ggsave(file='join_DimPlot_by4groups.jpg')
ggsave(file='join_DimPlot_by4groups.pdf')

library(dplyr)
#meta.data.table.sampled.byGroup <- meta.data.table %>% group_by(group) %>% sample_n(17165)

dat.CTRL_young <- subset(meta.data.table,group=="CTRL_young")[,c('UMAP_1','UMAP_2')]
dat.CTRL_old <- subset(meta.data.table,group=="CTRL_old")[,c('UMAP_1','UMAP_2')]
dat.CAD_without_T2D <- subset(meta.data.table,group=="CAD_without_T2D")[,c('UMAP_1','UMAP_2')]
dat.CAD_with_T2D <- subset(meta.data.table,group=="CAD_with_T2D")[,c('UMAP_1','UMAP_2')]

##### calculate density

get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

dat.CTRL_young$density <- get_density(dat.CTRL_young$UMAP_1, dat.CTRL_young$UMAP_2,n=100)
dat.CTRL_old$density <- get_density(dat.CTRL_old$UMAP_1, dat.CTRL_old$UMAP_2,n=100)
dat.CAD_without_T2D$density <- get_density(dat.CAD_without_T2D$UMAP_1, dat.CAD_without_T2D$UMAP_2,n=100)
dat.CAD_with_T2D$density <- get_density(dat.CAD_with_T2D$UMAP_1, dat.CAD_with_T2D$UMAP_2,n=100)

### contour map and gradiant color
library(viridis)
ggplot(dat.CTRL_young,aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(x=UMAP_1,y=UMAP_2, color = density)) +
 scale_color_viridis() + labs(x="UMAP 1",y='UMAP 2') + ggtitle('CTRL_young') + theme_bw()
ggsave(file="CTRL_young_densityPlot.jpg",width=4.75,height=3.68,units='in',dpi=800)

ggplot(dat.CTRL_old,aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(x=UMAP_1,y=UMAP_2, color = density)) +
 scale_color_viridis() + labs(x="UMAP 1",y='UMAP 2') + ggtitle('CTRL_old') + theme_bw()
ggsave(file="CTRL_old_densityPlot.jpg",width=4.75,height=3.68,units='in',dpi=800)

ggplot(dat.CAD_without_T2D,aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(x=UMAP_1,y=UMAP_2, color = density)) +
 scale_color_viridis() + labs(x="UMAP 1",y='UMAP 2') + ggtitle('CAD_without_T2D') + theme_bw()
ggsave(file="CAD_without_T2D_densityPlot.jpg",width=4.75,height=3.68,units='in',dpi=800)

ggplot(dat.CAD_with_T2D,aes(x=UMAP_1,y=UMAP_2)) + geom_point(aes(x=UMAP_1,y=UMAP_2, color = density)) +
 scale_color_viridis() + labs(x="UMAP 1",y='UMAP 2') + ggtitle('CAD_with_T2D') + theme_bw()
ggsave(file="CAD_with_T2D_densityPlot.jpg",width=4.75,height=3.68,units='in',dpi=800)

# proportional data prep ---------------
library(reshape2)

#sample lineage Level
obs.counts.sample.seurat_clusters <- as.data.frame(as.matrix(table(meta.data.table[,c('orig.ident','subcluster')])))
obs.counts.sample.seurat_clusters <- dcast(obs.counts.sample.seurat_clusters, orig.ident~subcluster)
write.table(obs.counts.sample.seurat_clusters,file='obs.counts.sample.seurat_clusters.tsv',row.names=F,col.names=T,quote=F,sep='\t')

obs.counts = as.matrix(read.table(file="obs.counts.sample.seurat_clusters.tsv", sep='\t',row.names = 1,header=T))
prop <- obs.counts/apply(obs.counts, 1, sum)*100
library(reshape2)
prop.melt.df <- melt(prop)
colnames(prop.melt.df) <- c('Sample','Cluster','Prop')
library(ggplot2)

library(plyr)
sampleNames <- c("C11","C14","C15","C16","C01","C02","C04","C05","C06","C07","C08","S02","S03","S05","S08",
"S11","S12","S13","S14","S16","S17","S18","S19","S20","S21","S22","S23","S24","S25","S27")

prop.melt.df$group <- mapvalues(as.character(prop.melt.df$Sample),from=sampleNames,
to=c("CTRL_old","CTRL_old","CTRL_old","CTRL_old","CTRL_old","CTRL_young","CTRL_young","CTRL_young","CTRL_old","CTRL_old","CTRL_old","CAD_without_T2D","CAD_without_T2D","CAD_without_T2D","CAD_with_T2D","CAD_with_T2D",
"CAD_without_T2D","CAD_without_T2D","CAD_with_T2D","CAD_with_T2D","CAD_without_T2D","CAD_without_T2D","CAD_with_T2D","CAD_with_T2D",
"CAD_with_T2D","CAD_without_T2D","CAD_with_T2D","CAD_with_T2D","CAD_with_T2D","CAD_with_T2D")
)

prop.melt.df$group <- factor(prop.melt.df$group,levels=c('CTRL_young','CTRL_old','CAD_without_T2D','CAD_with_T2D'))

library(dplyr)
library(plotrix) # std.error()
prop.melt.df_by_group <- prop.melt.df%>%group_by(Cluster,group)%>%summarise(Ave = mean(Prop),Ste = std.error(Prop))
#prop.melt.df_by_group$Cluster <- factor(prop.melt.df_by_subgroup$Cluster,levels=names(prop.diff))
# colored by subgroup
ggplot(prop.melt.df_by_group, aes(x=Cluster, y=Ave, fill=group)) +
geom_bar(stat = "identity", position=position_dodge(),alpha=0.8)  +
geom_errorbar(data=prop.melt.df_by_group,aes(ymin=Ave-Ste, ymax=Ave+Ste,color=group),
size=.8, width=.3, position=position_dodge(.9),alpha=1) +
#geom_point(data=prop.melt.df,aes(x=Cluster, y=Prop, color=group),position = position_jitter(width = 0.1),alpha=1,size=1) +
scale_fill_manual(values=c('#1f77b4','#ff7f0e','#2ca02c','#d62728')) +
scale_color_manual(values=c('#1f77b4','#ff7f0e','#2ca02c','#d62728')) +
theme_bw() + theme(legend.position = "top",legend.title=element_blank()) + 
labs(x='',y='Relative proportion (%)') +
 ylim(0,65) + scale_y_sqrt()
ggsave(file='obs.counts.sample.cellType.4groups.barplot.pdf')
ggsave(file='obs.counts.sample.cellType.4groups.barplot.ysqrt.pdf')


#-------- prop in each sample
library(ggplot2)
categary.fill.pallet  <- as.character(c(
  '#1f77b4','#ff7f0e','#2ca02c','#d62728','#7733B7',
  '#8c565B','#E236AF','#8C8C00','#700B35','#06BACE',
  "#193C3E", "#992756", "#DADD00", "#F76AA7", "#C4796A",
  "#246A73", "#FF6360", "#4DDD30", "#FFA449", "#74A6E8",
  "#955BC1","#FF0044","#F6757A","#265C42","#00D8FF",
  "#63C74D","#800026","#4567EA", "#E82C61","#5B74F4"))
sampleNames <- c("C11","C14","C15","C16","C01","C02","C04","C05","C06","C07","C08","S02","S03","S05","S08",
"S11","S12","S13","S14","S16","S17","S18","S19","S20","S21","S22","S23","S24","S25","S27")

prop.melt.df$Sample <- factor(prop.melt.df$Sample,levels=c('C02','C04','C05',
'C01','C06','C07','C08','C11','C14','C15','C16',
'S02','S03','S05','S12','S13','S17','S18','S22',
'S08','S11','S14','S16','S19','S20','S21','S23','S24','S25','S27'))

prop.melt.df$Prop <- prop.melt.df$Prop*100
ggplot(prop.melt.df, aes(x=Sample, y=Prop, fill=Cluster)) + geom_bar(stat = "identity",position = "fill") +
theme_bw() +labs(x='',y='Relative proportion (%)') + scale_fill_manual(values=categary.fill.pallet)
ggsave(file='cellType_prop_eachSample.pdf')


df <- subset(prop.melt.df,Cluster=="Ad2")
res.aov <- aov(Prop ~ group, data = df)
summary(res.aov)
            Df   Sum Sq Mean Sq F value  Pr(>F)
group        3  8288502 2762834   4.976 0.00735 **
Residuals   26 14436791  555261
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


library(multcomp)
summary(glht(res.aov, linfct = mcp(group = "Tukey")))

Fit: aov(formula = Prop ~ subgroup, data = df)
Linear Hypotheses:
                                    Estimate Std. Error t value Pr(>|t|)
CTRL_old - CTRL_young == 0          -835.281    504.475  -1.656   0.3604
CAD_without_T2D - CTRL_young == 0    593.767    504.475   1.177   0.6407
CAD_with_T2D - CTRL_young == 0         7.424    485.351   0.015   1.0000
CAD_without_T2D - CTRL_old == 0     1429.049    372.579   3.836   0.0036 **
CAD_with_T2D - CTRL_old == 0         842.705    346.245   2.434   0.0928 .
CAD_with_T2D - CAD_without_T2D == 0 -586.344    346.245  -1.693   0.3412

