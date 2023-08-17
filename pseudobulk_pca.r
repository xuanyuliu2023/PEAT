library(Seurat)
library(DESeq2)
load('/work/project/xuanyu/project/PEAT_V2/seurat/integrated/meta.data.table.RData')
load('/work/project/xuanyu/project/PEAT_V2/seurat/integrated/seuObj.integrated.ALL.RData')

# 体量的稀疏矩阵转换为稠密矩阵，使用C++函数 as_matrix
library(Rcpp)
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){
  int k = z.size() ;
  IntegerMatrix  mat(nrows, ncols);
  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }
  return mat;
}
' )

as_matrix <- function(mat){
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

counts <- as.data.frame(as_matrix(seuObj.integrated@assays$RNA@counts))


Sample1.cell <- row.names(subset(meta.data.table,orig.ident=='C01' ))
Sample2.cell <- row.names(subset(meta.data.table,orig.ident=='C06' ))
Sample3.cell <- row.names(subset(meta.data.table,orig.ident=='C07' ))
Sample4.cell <- row.names(subset(meta.data.table,orig.ident=='C08' ))
Sample5.cell <- row.names(subset(meta.data.table,orig.ident=='C11' ))
Sample6.cell <- row.names(subset(meta.data.table,orig.ident=='C14' ))
Sample7.cell <- row.names(subset(meta.data.table,orig.ident=='C15' ))
Sample8.cell <- row.names(subset(meta.data.table,orig.ident=='C16' ))
Sample9.cell <- row.names(subset(meta.data.table,orig.ident=='S02' ))
Sample10.cell <- row.names(subset(meta.data.table,orig.ident=='S03'))
Sample11.cell <- row.names(subset(meta.data.table,orig.ident=='S05'))
Sample12.cell <- row.names(subset(meta.data.table,orig.ident=='S12'))
Sample13.cell <- row.names(subset(meta.data.table,orig.ident=='S13'))
Sample14.cell <- row.names(subset(meta.data.table,orig.ident=='S17'))
Sample15.cell <- row.names(subset(meta.data.table,orig.ident=='S18'))
Sample16.cell <- row.names(subset(meta.data.table,orig.ident=='S22'))
Sample17.cell <- row.names(subset(meta.data.table,orig.ident=='S08'))
Sample18.cell <- row.names(subset(meta.data.table,orig.ident=='S11'))
Sample19.cell <- row.names(subset(meta.data.table,orig.ident=='S14'))
Sample20.cell <- row.names(subset(meta.data.table,orig.ident=='S16'))
Sample21.cell <- row.names(subset(meta.data.table,orig.ident=='S19'))
Sample22.cell <- row.names(subset(meta.data.table,orig.ident=='S20'))
Sample23.cell <- row.names(subset(meta.data.table,orig.ident=='S21'))
Sample24.cell <- row.names(subset(meta.data.table,orig.ident=='S25'))


Sample1.counts <- counts[,Sample1.cell]
Sample2.counts <- counts[,Sample2.cell]
Sample3.counts <- counts[,Sample3.cell]
Sample4.counts <- counts[,Sample4.cell]
Sample5.counts <- counts[,Sample5.cell]
Sample6.counts <- counts[,Sample6.cell]
Sample7.counts <- counts[,Sample7.cell]
Sample8.counts <- counts[,Sample8.cell]
Sample9.counts <- counts[,Sample9.cell]
Sample10.counts <- counts[,Sample10.cell]
Sample11.counts <- counts[,Sample11.cell]
Sample12.counts <- counts[,Sample12.cell]
Sample13.counts <- counts[,Sample13.cell]
Sample14.counts <- counts[,Sample14.cell]
Sample15.counts <- counts[,Sample15.cell]
Sample16.counts <- counts[,Sample16.cell]
Sample17.counts <- counts[,Sample17.cell]
Sample18.counts <- counts[,Sample18.cell]
Sample19.counts <- counts[,Sample19.cell]
Sample20.counts <- counts[,Sample20.cell]
Sample21.counts <- counts[,Sample21.cell]
Sample22.counts <- counts[,Sample22.cell]
Sample23.counts <- counts[,Sample23.cell]
Sample24.counts <- counts[,Sample24.cell]

library(Matrix)
Sample1.data <- rowSums(Sample1.counts)
Sample2.data <- rowSums(Sample2.counts)
Sample3.data <- rowSums(Sample3.counts)
Sample4.data <- rowSums(Sample4.counts)
Sample5.data <- rowSums(Sample5.counts)
Sample6.data <- rowSums(Sample6.counts)
Sample7.data <- rowSums(Sample7.counts)
Sample8.data <- rowSums(Sample8.counts)
Sample9.data <- rowSums(Sample9.counts)
Sample10.data <- rowSums(Sample10.counts)
Sample11.data <- rowSums(Sample11.counts)
Sample12.data <- rowSums(Sample12.counts)
Sample13.data <- rowSums(Sample13.counts)
Sample14.data <- rowSums(Sample14.counts)
Sample15.data <- rowSums(Sample15.counts)
Sample16.data <- rowSums(Sample16.counts)
Sample17.data <- rowSums(Sample17.counts)
Sample18.data <- rowSums(Sample18.counts)
Sample19.data <- rowSums(Sample19.counts)
Sample20.data <- rowSums(Sample20.counts)
Sample21.data <- rowSums(Sample21.counts)
Sample22.data <- rowSums(Sample22.counts)
Sample23.data <- rowSums(Sample23.counts)
Sample24.data <- rowSums(Sample24.counts)


df <- data.frame(C01=Sample1.data,C06=Sample2.data,C07=Sample3.data,C08=Sample4.data,
                 C11=Sample5.data,C14=Sample6.data,C15=Sample7.data,C16=Sample8.data,
                S02=Sample9.data,S03=Sample10.data,S05=Sample11.data,S12=Sample12.data,
                S13=Sample13.data,S17=Sample14.data,S18=Sample15.data,S22=Sample16.data,
                S08=Sample17.data,S11=Sample18.data,S14=Sample19.data,S16=Sample20.data,
                S19=Sample21.data,S20=Sample22.data,S21=Sample23.data,S25=Sample24.data)

library(plyr)
subgroup <- mapvalues(colnames(df),from= as.character(meta.data.table$orig.ident),to= as.character(meta.data.table$group))
sex <- mapvalues(colnames(df),from= as.character(meta.data.table$orig.ident),to= as.character(meta.data.table$sex))
coldata <- data.frame(condition =subgroup, row.names=colnames(df))
coldata$sex <- sex
library("BiocParallel")
register(MulticoreParam(10))


coldata_subset = coldata
sampleNames =rownames(coldata_subset)
geneType <- read.table(file='/work/project/xuanyu/resource/10XGenomics/refdata-gex-GRCh38-2020-A/genes.gtf.geneType.tsv',sep='\t',header=F)
protein_coding_genes <- subset(geneType, V2 == 'protein_coding')$V1
df <- df[row.names(df) %in% protein_coding_genes,]
df_subset = df[,sampleNames]

#By adding variables to the design, one can control for additional variation in the counts. For example, if the condition samples are balanced across experimental batches, by including the batch factor to the design, one can increase the sensitivity for finding differences due to condition.

dds <- DESeqDataSetFromMatrix(countData = df_subset,
                              colData = coldata_subset,
                              design= ~ sex + condition)
save(coldata_subset,file='coldata_subset.Rdata')
#Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set the reference level
dds$condition <- factor(dds$condition, levels = c('CTRL','CAD_without_T2D','CAD_with_T2D'))

dds <- DESeq(dds)
plotCounts(dds, gene='ADIPOQ', intgroup="condition")

# Extracting transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
save(vsd,file="vsd.Rdata")















# gene expression plot
library(reshape2)
library(ggplot2)
library(DESeq2)
library(ggrepel)

load("vsd.Rdata")
load('coldata_subset.Rdata')
normalizedCountMatrix <- assay(vsd)


df_t <- t(normalizedCountMatrix)
#var.genes <- seuObj.integrated@assays$RNA@var.features
#df_t_var.genes <- df_t[,var.genes]

df.pca <- prcomp(df_t,center = TRUE,scale. = TRUE)
pca.score <- cbind(as.data.frame(df.pca$x[,c('PC1','PC2','PC3','PC4')]),sample=rownames(df_t))

library(plyr)
pca.score$group <- mapvalues(rownames(pca.score),from=as.character(meta.data.table$orig.ident),to=as.character(meta.data.table$group))
pca.score$sex <- mapvalues(rownames(pca.score),from=as.character(meta.data.table$orig.ident),to=as.character(meta.data.table$sex))
pca.score$sample <- row.names(pca.score)
pc1PV <- round(summary(df.pca)$importance[2,1]*100,2)
pc2PV <- round(summary(df.pca)$importance[2,2]*100,2)
pc3PV <- round(summary(df.pca)$importance[2,3]*100,2)
pc4PV <- round(summary(df.pca)$importance[2,4]*100,2)
xLAB <- paste("PC1 (",pc1PV,"% of Variance Explained)",sep='')
yLAB <- paste("PC2 (",pc2PV,"% of Variance Explained)",sep='')
zLAB <- paste("PC3 (",pc3PV,"% of Variance Explained)",sep='')

library(ggrepel)
library(ggplot2)
my.plot <- ggplot(pca.score,aes(PC1,PC2,color=group,label=sample))
my.plot + geom_point(size=3,alpha=1) + geom_hline(yintercept = 0,linetype="dotdash",alpha=0.5) +
 geom_vline(xintercept = 0,linetype="dotdash",alpha=0.3) +
 geom_text_repel(size=3) +
 labs(x=xLAB,y=yLAB)+ theme_bw() +
 stat_ellipse(type ="norm",linetype=2,level=0.4) +
 theme(legend.title=element_blank(), legend.position='top') + scale_color_manual(values = rev(c('#2ca02c','#ff7f0e','#d62728')))
ggsave('allCells_pca_PC1_PC2.pdf')

my.plot <- ggplot(pca.score,aes(PC1,PC2,color=sex,label=sex))
my.plot + geom_point(size=3,alpha=1) + geom_hline(yintercept = 0,linetype="dotdash",alpha=0.5) +
 geom_vline(xintercept = 0,linetype="dotdash",alpha=0.3) +
 geom_text_repel(size=3) +
 labs(x=xLAB,y=yLAB)+ theme_bw() +
 theme(legend.title=element_blank(), legend.position='top') + scale_color_manual(values = rev(c('#1f77b4','#d62728')))
ggsave('allCells_pca_PC1_PC2_sex.pdf')

library(ggrepel)
library(ggplot2)
my.plot <- ggplot(pca.score,aes(PC3,PC4,color=group,label=sample))
my.plot + geom_point(size=3,alpha=1) + geom_hline(yintercept = 0,linetype="dotdash",alpha=0.5) +
 geom_vline(xintercept = 0,linetype="dotdash",alpha=0.3) +
 geom_text_repel(size=3) +
 labs(x=xLAB,y=yLAB)+ theme_bw() +
 theme(legend.title=element_blank(), legend.position='top') + scale_color_manual(values = rev(c('#2ca02c','#ff7f0e','#d62728')))
ggsave('allCells_pca_PC3_PC4.pdf')

library(plotly)

text.p = ~paste('<br>sample:', sample, 'subgroup: ',subgroup,'<br>sex:', sex)
fig <- plot_ly(pca.score, x = ~PC1, y = ~PC2, z = ~PC3, color = ~subgroup, colors =  rev(c('#2ca02c','#ff7f0e','#d62728')),text=text.p)

htmlwidgets::saveWidget(fig, "pca_3D.html")
