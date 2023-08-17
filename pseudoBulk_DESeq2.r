library(Seurat)
library(DESeq2)
load('/work/project/xuanyu/project/PEAT_V2/seurat/integrated/meta.data.table.RData')
load('/work/project/xuanyu/project/PEAT_V2/seurat/integrated/seuObj.integrated.ALL.RData')

celltype = "BC"
CTRL_group = 'CAD_without_T2D'
CASE_group = 'CAD_with_T2D'



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



Sample1.cell <- row.names(subset(meta.data.table,orig.ident=='C01' & cellType==celltype))
Sample2.cell <- row.names(subset(meta.data.table,orig.ident=='C06' & cellType==celltype))
Sample3.cell <- row.names(subset(meta.data.table,orig.ident=='C07' & cellType==celltype))
Sample4.cell <- row.names(subset(meta.data.table,orig.ident=='C08' & cellType==celltype))
Sample5.cell <- row.names(subset(meta.data.table,orig.ident=='C11' & cellType==celltype))
Sample6.cell <- row.names(subset(meta.data.table,orig.ident=='C14' & cellType==celltype))
Sample7.cell <- row.names(subset(meta.data.table,orig.ident=='C15' & cellType==celltype))
Sample8.cell <- row.names(subset(meta.data.table,orig.ident=='C16' & cellType==celltype))
Sample9.cell <- row.names(subset(meta.data.table,orig.ident=='S02' & cellType==celltype))
Sample10.cell <- row.names(subset(meta.data.table,orig.ident=='S03' & cellType==celltype))
Sample11.cell <- row.names(subset(meta.data.table,orig.ident=='S05' & cellType==celltype))
Sample12.cell <- row.names(subset(meta.data.table,orig.ident=='S12' & cellType==celltype))
Sample13.cell <- row.names(subset(meta.data.table,orig.ident=='S13' & cellType==celltype))
Sample14.cell <- row.names(subset(meta.data.table,orig.ident=='S17' & cellType==celltype))
Sample15.cell <- row.names(subset(meta.data.table,orig.ident=='S18' & cellType==celltype))
Sample16.cell <- row.names(subset(meta.data.table,orig.ident=='S22' & cellType==celltype))
Sample17.cell <- row.names(subset(meta.data.table,orig.ident=='S08' & cellType==celltype))
Sample18.cell <- row.names(subset(meta.data.table,orig.ident=='S11' & cellType==celltype))
Sample19.cell <- row.names(subset(meta.data.table,orig.ident=='S14' & cellType==celltype))
Sample20.cell <- row.names(subset(meta.data.table,orig.ident=='S16' & cellType==celltype))
Sample21.cell <- row.names(subset(meta.data.table,orig.ident=='S19' & cellType==celltype))
Sample22.cell <- row.names(subset(meta.data.table,orig.ident=='S20' & cellType==celltype))
Sample23.cell <- row.names(subset(meta.data.table,orig.ident=='S21' & cellType==celltype))
Sample24.cell <- row.names(subset(meta.data.table,orig.ident=='S25' & cellType==celltype))


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




coldata_subset = subset(coldata,condition %in% c(CTRL_group,CASE_group))
sampleNames =rownames(coldata_subset)
geneType <- read.table(file='/work/project/xuanyu/resource/10XGenomics/refdata-gex-GRCh38-2020-A/genes.gtf.geneType.tsv',sep='\t',header=F)
protein_coding_genes <- subset(geneType, V2 == 'protein_coding')$V1
df <- df[row.names(df) %in% protein_coding_genes,]
df_subset = df[,sampleNames]

#By adding variables to the design, one can control for additional variation in the counts. For example, if the condition samples are balanced across experimental batches, by including the batch factor to the design, one can increase the sensitivity for finding differences due to condition.

dds <- DESeqDataSetFromMatrix(countData = df_subset,
                              colData = coldata_subset,
                              design= ~  sex + condition)
save(coldata_subset,file='coldata_subset.Rdata')

#Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set the reference level
dds$condition <- factor(dds$condition, levels = c(CTRL_group,CASE_group))

dds <- DESeq(dds)
res <- results(dds, name=paste0("condition_",CASE_group,"_vs_",CTRL_group), alpha=0.05)
#res <- lfcShrink(dds, coef="condition_CASE_vs_CTRL", type="ashr")

summary(res)
sum(res$padj < 0.05, na.rm=TRUE)
res.df <- data.frame(res)
res.df$gene <- row.names(res.df)

res.df <- res.df[order(res.df$pvalue),]
res.sig.df <- subset(res.df,padj < 0.05 & abs(log2FoldChange) >1)
dim(res.sig.df)
write.table(res.df,file=paste0("condition_",CASE_group,"_vs_",CTRL_group,".pseudobulk_results.tsv"),sep= '\t',row.names=F,col.names=T,quote=F)
write.table(res.sig.df,file=paste0("condition_",CASE_group,"_vs_",CTRL_group,".pseudobulk_results.padj0.05.tsv"),sep= '\t',row.names=F,col.names=T,quote=F)

save(dds,res,file=paste0(celltype,'.res.Rdata'))

#plotCounts(dds, gene='MRC1', intgroup="condition")

# Extracting transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
save(vsd,file=paste0(celltype,"_vsd_.Rdata"))

