library(CellChat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(igraph)
# Create a directory to save figures
data.dir <- './comparison'
dir.create(data.dir)
setwd(data.dir)
#Load CellChat object of each dataset and then merge together
load('../cellchat.CTRL.RData')
load('../cellchat.CAD_without_T2D.RData')
load('../cellchat.CAD_with_T2D.RData')

# all three samples object
object.list <- list(CTRL = cellchat.CTRL, CAD_without_T2D = cellchat.CAD_without_T2D, CAD_with_T2D = cellchat.CAD_with_T2D)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))



# two samples object
object.list.CTRL.CAD_without_T2D <- list(CTRL = cellchat.CTRL, CAD_without_T2D = cellchat.CAD_without_T2D)
cellchat.CTRL.CAD_without_T2D <- mergeCellChat(object.list.CTRL.CAD_without_T2D, add.names = names(object.list.CTRL.CAD_without_T2D))

#save(cellchat,file='cellchat.merged.Rdata')

#Part I: Predict general principles of cell-cell communication
# Compare the total number of interactions and interaction strength
group.colors <- c('#2b9f2b','#ff7f0e','#d62526')
gg1 <- compareInteractions(cellchat, show.legend = F, color.use = group.colors, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F,color.use = group.colors, group = c(1,2,3), measure = "weight")
x11()
gg1 + gg2
dev.copy2pdf(file='totalNumberofInteractionsAndInteractionStrength.pdf')


#Compare the number of interactions and interaction strength among different cell populations
#x11()
#par(mfrow = c(1,2), xpd=TRUE)
#netVisual_diffInteraction(cellchat, weight.scale = T,remove.isolate=F,label.edge= T)
#netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",remove.isolate=F,label.edge= T)
#dev.copy2pdf(file='DifferentialNumberofInteractionsAndInteractionStrength.circlePlot.pdf')


# heatmap
# The top colored bar plot represents the sum of column of values displayed in the heatmap (incoming signaling).
# The right colored bar plot represents the sum of row of values (outgoing signaling).
# In the colorbar, red or blue represents increased or decreased signaling in the second dataset compared to the first one.
x11()
gg1 <- netVisual_heatmap(cellchat,comparison = c(1, 2),remove.isolate=F,cluster.rows=F,cluster.cols=F)
gg1
dev.copy2pdf(file='CAD_without_T2D_vs_CTRL_totalNumberofInteractions.heatmap.pdf')
x11()
gg1 <- netVisual_heatmap(cellchat,comparison = c(1, 3),remove.isolate=F,cluster.rows=F,cluster.cols=F)
gg1
dev.copy2pdf(file='CAD_with_T2D_vs_CTRL_totalNumberofInteractions.heatmap.pdf')
x11()
gg1 <- netVisual_heatmap(cellchat,comparison = c(2, 3),remove.isolate=F,cluster.rows=F,cluster.cols=F)
gg1
dev.copy2pdf(file='CAD_with_T2D_vs_CAD_without_T2D_totalNumberofInteractions.heatmap.pdf')


#> Do heatmap showing InteractionStrength based on a merged object
x11()
gg1 <- netVisual_heatmap(cellchat,comparison = c(1, 2),remove.isolate=F,cluster.rows=F,cluster.cols=F,measure = "weight")
gg1
dev.copy2pdf(file='CAD_without_T2D_vs_CTRL_InteractionStrength.heatmap.pdf')
x11()
gg1 <- netVisual_heatmap(cellchat,comparison = c(1, 3),remove.isolate=F,cluster.rows=F,cluster.cols=F,measure = "weight")
gg1
dev.copy2pdf(file='CAD_with_T2D_vs_CTRL_InteractionStrength.heatmap.pdf')
x11()
gg1 <- netVisual_heatmap(cellchat,comparison = c(2, 3),remove.isolate=F,cluster.rows=F,cluster.cols=F,measure = "weight")
gg1
dev.copy2pdf(file='CAD_with_T2D_vs_CAD_without_T2D_InteractionStrength.heatmap.pdf')



weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
x11()
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, vertex.weight=as.numeric(table(object.list[[i]]@idents)),remove.isolate=F, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2],
   edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.copy2pdf(file='totalNumberofInteractions.circlePlot.pdf')


#Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) + theme_bw()
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways

x11()
patchwork::wrap_plots(plots = gg)
dev.copy2pdf(file='interactionStrength2D.scatter.pdf')

#------------------DifferentialInteractionStrength2D.scatter

NL_centrality_df <- gg[[1]]$data
LS_centrality_df <- gg[[2]]$data
x_cord <- LS_centrality_df$x - NL_centrality_df$x
y_cord <- LS_centrality_df$y - NL_centrality_df$y
diff_centrality_df <- data.frame(x=x_cord,y=y_cord,labels=NL_centrality_df$labels)
library(ggplot2)
library(ggrepel)
color.p <- c('#e41a1c','#377eb8','#4daf4a','#984ea3',
'#f29403','#f781bf','#bc9dcc','#a65628',
'#54b0e4','#222f75','#1b9e77','#b2df8a',
'#e3be00','#fb9a99','#e7298a')
ggplot(diff_centrality_df,aes(x=x,y=y,color=labels,label=labels)) + geom_point(size=5,alpha=1) + geom_hline(yintercept = 0,linetype="dotdash",alpha=0.5) +
 geom_vline(xintercept = 0,linetype="dotdash",alpha=0.3) +
 geom_text_repel(size=3) +
 labs(x="Differential outgoing interaction strength",y="Differential incoming interaction strength")+ theme_bw() +
 theme(legend.title=element_blank(), legend.justification=c(-0.1,1.1), legend.position="none") + scale_color_manual(values=color.p)
ggsave(file="CAD_without_T2D_vs_CTRL_DifferentialInteractionStrength2D.scatter.pdf")


NL_centrality_df <- gg[[1]]$data
LS_centrality_df <- gg[[3]]$data
x_cord <- LS_centrality_df$x - NL_centrality_df$x
y_cord <- LS_centrality_df$y - NL_centrality_df$y
diff_centrality_df <- data.frame(x=x_cord,y=y_cord,labels=NL_centrality_df$labels)
library(ggplot2)
library(ggrepel)
color.p <- c('#e41a1c','#377eb8','#4daf4a','#984ea3',
'#f29403','#f781bf','#bc9dcc','#a65628',
'#54b0e4','#222f75','#1b9e77','#b2df8a',
'#e3be00','#fb9a99','#e7298a')
ggplot(diff_centrality_df,aes(x=x,y=y,color=labels,label=labels)) + geom_point(size=5,alpha=1) + geom_hline(yintercept = 0,linetype="dotdash",alpha=0.5) +
 geom_vline(xintercept = 0,linetype="dotdash",alpha=0.3) +
 geom_text_repel(size=3) +
 labs(x="Differential outgoing interaction strength",y="Differential incoming interaction strength")+ theme_bw() +
 theme(legend.title=element_blank(), legend.justification=c(-0.1,1.1), legend.position="none") + scale_color_manual(values=color.p)
ggsave(file="CAD_with_T2D_vs_CTRL_DifferentialInteractionStrength2D.scatter.pdf")


NL_centrality_df <- gg[[2]]$data
LS_centrality_df <- gg[[3]]$data
x_cord <- LS_centrality_df$x - NL_centrality_df$x
y_cord <- LS_centrality_df$y - NL_centrality_df$y
diff_centrality_df <- data.frame(x=x_cord,y=y_cord,labels=NL_centrality_df$labels)
library(ggplot2)
library(ggrepel)
color.p <- c('#e41a1c','#377eb8','#4daf4a','#984ea3',
'#f29403','#f781bf','#bc9dcc','#a65628',
'#54b0e4','#222f75','#1b9e77','#b2df8a',
'#e3be00','#fb9a99','#e7298a')
ggplot(diff_centrality_df,aes(x=x,y=y,color=labels,label=labels)) + geom_point(size=5,alpha=1) + geom_hline(yintercept = 0,linetype="dotdash",alpha=0.5) +
 geom_vline(xintercept = 0,linetype="dotdash",alpha=0.3) +
 geom_text_repel(size=3) +
 labs(x="Differential outgoing interaction strength",y="Differential incoming interaction strength")+ theme_bw() +
 theme(legend.title=element_blank(), legend.justification=c(-0.1,1.1), legend.position="none") + scale_color_manual(values=color.p)
ggsave(file="CAD_with_T2D_vs_CAD_without_T2D_DifferentialInteractionStrength2D.scatter.pdf")


# Identify signaling changes associated with one cell group
# Please note that the title of each plot maybe wrong, e.g. CTRL_old_ vs CAD_without_T2D
x11()
cellType = "Adipocyte"
gg1 <- netAnalysis_signalingChanges_scatter(object.list,comparison = c(1, 2), idents.use = cellType)
gg2 <- netAnalysis_signalingChanges_scatter(object.list,comparison = c(1, 3), idents.use = cellType)
gg3 <- netAnalysis_signalingChanges_scatter(object.list,comparison = c(2, 3), idents.use = cellType)
wrap_plots(gg1,gg2,gg3)
dev.copy2pdf(file=paste0(cellType,'_Differential_interactionStrength2D.scatter.pdf'))


#Part II: Identify the conserved and context-specific signaling pathways
#Identify signaling networks with larger (or less) difference as well as
#signaling groups based on their functional/structure similarity

# Identify signaling groups based on their functional similarity
# Functional similarity: High degree of functional similarity indicates major senders and receivers are similar,
# and it can be interpreted as the two signaling pathways or two ligand-receptor pairs exhibit
# similar and/or redundant roles
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional",umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
x11()
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5) + theme_bw()
dev.copy2pdf(file="signalingGroupFunctionalSimilariy")

# A structural similarity was used to compare their signaling network structure, without considering the similarity of senders and receivers
# slow!
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural",umap.method = "uwot")
cellchat <- netClustering(cellchat, type = "structural")

# Visualization in 2D-space
x11()
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5)
dev.copy2pdf(file='signalingGroupStructuralSimilariy.pdf')

# Compute and visualize the pathway distance in the learned joint manifold
#Larger distance implies larger difference of the communication networks between two datasets in terms
#of either functional or structure similarity

x11()
rankSimilarity(cellchat, type = "functional",comparison2 = c(1, 2))
dev.copy2pdf(file="CAD_without_T2D_vs_CTRL.pathwayDistance.pdf")

x11()
rankSimilarity(cellchat, type = "functional",comparison2 = c(1, 3))
dev.copy2pdf(file="CAD_with_T2D_vs_CTRL.pathwayDistance.pdf")

x11()
rankSimilarity(cellchat, type = "functional",comparison2 = c(2, 3))
dev.copy2pdf(file="CAD_with_T2D_vs_CAD_without_T2D.pathwayDistance.pdf")

#Identify and visualize the conserved and context-specific signaling pathways
#By comparing the information flow/interaction strengh of each signaling pathway,
# we can identify signalin137sampsJan6g pathways, (i) turn off, (ii) decrease, (iii) turn on or (iv) increase,
#by change their information flow at one condition as compared to another condition.
# Compare the overall information flow of each signaling pathway
# Information flow for each signaling pathway, which is defined by the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network).
group.colors <- c('#2b9f2b','#ff7f0e','#d62526')
x11()
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, comparison = c(1, 2),color.use =c('#2b9f2b','#ff7f0e'))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, comparison = c(1, 2),color.use =c('#2b9f2b','#ff7f0e'))
gg1 + gg2
dev.copy2pdf(file='CAD_without_T2D_vs_CTRL.InformationFlow.pdf')

x11()
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, comparison = c(1, 3),color.use =c('#2b9f2b','#d62526'))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, comparison = c(1, 3),color.use =c('#2b9f2b','#d62526'))
gg1 + gg2
dev.copy2pdf(file='CAD_with_T2D_vs_CTRL_old.InformationFlow.pdf')

x11()
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, comparison = c(2, 3),color.use =c('#ff7f0e','#d62526'))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, comparison = c(2, 3),color.use =c('#ff7f0e','#d62526'))
gg1 + gg2
dev.copy2pdf(file='CAD_with_T2D_vs_CAD_without_T2D.InformationFlow.pdf')


# all 
x11()
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, comparison = c(1,2,3),color.use =c('#2b9f2b','#ff7f0e','#d62526'))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, comparison = c(1,2,3),color.use =c('#2b9f2b','#ff7f0e','#d62526'))
gg1 + gg2
dev.copy2pdf(file='ALL_groups.InformationFlow.pdf')


# Compare outgoing (or incoming) signaling associated with each cell population
library(ComplexHeatmap)
object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]])
object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]])
object.list[[3]] <- netAnalysis_computeCentrality(object.list[[3]])


# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, 
title = names(object.list)[1], width = 9, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, 
title = names(object.list)[2], width = 9, height = 20)
x11()
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='CAD_without_T2D_vs_CTRL_outgoingSignalingPathways.pdf')

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, 
title = names(object.list)[1], width = 9, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, 
title = names(object.list)[2], width = 9, height = 20, color.heatmap = "GnBu")
x11()
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='CAD_without_T2D_vs_CTRL_IncomingSignalingPathways.pdf')

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", signaling = pathway.union, 
title = names(object.list)[1], width = 9, height = 20, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", signaling = pathway.union, 
title = names(object.list)[2], width = 9, height = 20, color.heatmap = "OrRd")
x11()
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='CAD_without_T2D_vs_CTRL_OverallSignalingPathways.pdf')
# ---------
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[3]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, 
title = names(object.list)[1], width = 9, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "outgoing", signaling = pathway.union, 
title = names(object.list)[3], width = 9, height = 20)
x11()
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='CAD_with_T2D_vs_CTRL_outgoingSignalingPathways.pdf')

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, 
title = names(object.list)[1], width = 9, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "incoming", signaling = pathway.union, 
title = names(object.list)[3], width = 9, height = 20, color.heatmap = "GnBu")
x11()
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='CAD_with_T2D_vs_CTRL_IncomingSignalingPathways.pdf')

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", signaling = pathway.union, 
title = names(object.list)[1], width = 9, height = 20, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "all", signaling = pathway.union, 
title = names(object.list)[3], width = 9, height = 20, color.heatmap = "OrRd")
x11()
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='CAD_with_T2D_vs_CTRL_OverallSignalingPathways.pdf')
#------
pathway.union <- union(object.list[[2]]@netP$pathways, object.list[[3]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, 
title = names(object.list)[2], width = 9, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "outgoing", signaling = pathway.union, 
title = names(object.list)[3], width = 9, height = 20)
x11()
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='CAD_with_T2D_vs_CAD_without_T2D_outgoingSignalingPathways.pdf')

ht1 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, 
title = names(object.list)[2], width = 9, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "incoming", signaling = pathway.union, 
title = names(object.list)[3], width = 9, height = 20, color.heatmap = "GnBu")
x11()
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='CAD_with_T2D_vs_CAD_without_T2D_IncomingSignalingPathways.pdf')

ht1 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", signaling = pathway.union, 
title = names(object.list)[2], width = 9, height = 20, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "all", signaling = pathway.union, 
title = names(object.list)[3], width = 9, height = 20, color.heatmap = "OrRd")
x11()
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='CAD_with_T2D_vs_CAD_without_T2D_OverallSignalingPathways.pdf')

# --- plot all together
pathway.union <- union(union(object.list[[1]]@netP$pathways,object.list[[2]]@netP$pathways),object.list[[3]]@netP$pathways)

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = pathway.union, 
title = names(object.list)[1], width = 9, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = pathway.union, 
title = names(object.list)[2], width = 9, height = 20)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "outgoing", signaling = pathway.union, 
title = names(object.list)[3], width = 9, height = 20)

x11()
draw(ht1 + ht2 + ht3 , ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='ALL_outgoingSignalingPathways.pdf')

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, 
title = names(object.list)[1], width = 9, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, 
title = names(object.list)[2], width = 9, height = 20, color.heatmap = "GnBu")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "incoming", signaling = pathway.union, 
title = names(object.list)[3], width = 9, height = 20, color.heatmap = "GnBu")

x11()
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='ALL_IncomingSignalingPathways.pdf')

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", signaling = pathway.union, 
title = names(object.list)[1], width = 9, height = 20, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", signaling = pathway.union, 
title = names(object.list)[2], width = 9, height = 20, color.heatmap = "OrRd")
ht3 = netAnalysis_signalingRole_heatmap(object.list[[3]], pattern = "all", signaling = pathway.union, 
title = names(object.list)[3], width = 9, height = 20, color.heatmap = "OrRd")

x11()
draw(ht1 + ht2 + ht3, ht_gap = unit(0.5, "cm"))
dev.copy2pdf(file='ALL_OverallSignalingPathways.pdf')

#-----save
save(cellchat,object.list,file='cellchat_object.list.Rdata')



#Part III: Identify the upgulated and down-regulated signaling ligand-receptor pairs
#------- Method1
# Identify dysfunctional signaling by comparing the communication probabities
source.cellType <- c('Adipocyte')
target.cellType <- c('Macrophage','ASPC','Adipocyte')
NL_groupName <- c('CTRL')
LS_groupName <- c('CAD_without_T2D')

source.cellType <- c('Adipocyte')
target.cellType <- c('Macrophage','ASPC','Adipocyte')
NL_groupName <- c('CTRL')
LS_groupName <- c('CAD_with_T2D')

source.cellType <- c('Adipocyte')
target.cellType <- c('Macrophage','ASPC','Adipocyte')
NL_groupName <- c('CAD_without_T2D')
LS_groupName <- c('CAD_with_T2D')


source.cellType <- c('ASPC')
target.cellType <- c('Macrophage','ASPC','Adipocyte')
NL_groupName <- c('CTRL')
LS_groupName <- c('CAD_without_T2D')

source.cellType <- c('ASPC')
target.cellType <- c('Macrophage','ASPC','Adipocyte')
NL_groupName <- c('CTRL')
LS_groupName <- c('CAD_with_T2D')

source.cellType <- c('ASPC')
target.cellType <- c('Macrophage','ASPC','Adipocyte')
NL_groupName <- c('CAD_without_T2D')
LS_groupName <- c('CAD_with_T2D')


source.cellType <- c('Macrophage')
target.cellType <- c('Macrophage','ASPC','Adipocyte')
NL_groupName <- c('CTRL')
LS_groupName <- c('CAD_without_T2D')

source.cellType <- c('Macrophage')
target.cellType <- c('Macrophage','ASPC','Adipocyte')
NL_groupName <- c('CTRL')
LS_groupName <- c('CAD_with_T2D')

source.cellType <- c('Macrophage')
target.cellType <- c('Macrophage','ASPC','Adipocyte')
NL_groupName <- c('CAD_without_T2D')
LS_groupName <- c('CAD_with_T2D')


dataSetVec <- c(1,2,3)
names(dataSetVec) <- c('CTRL','CAD_without_T2D','CAD_with_T2D')
gg1 <- netVisual_bubble(cellchat, sources.use = source.cellType, targets.use = target.cellType,  
 comparison = c(dataSetVec[NL_groupName],dataSetVec[LS_groupName]), max.dataset = dataSetVec[LS_groupName],
  title.name = paste0("Increased signaling in ",LS_groupName,' vs ',NL_groupName), 
 angle.x = 45, remove.isolate = T,color.heatmap='Spectral')
 
#gg1$data <- gg1$data[gg1$data$interaction_name %in%  net.up.CTRL.CAD_without_T2D$interaction_name,]

gg2 <- netVisual_bubble(cellchat, sources.use = source.cellType, targets.use = target.cellType,
comparison = c(dataSetVec[NL_groupName],dataSetVec[LS_groupName]), max.dataset = dataSetVec[NL_groupName], 
title.name = paste0("Decreased signaling in ",LS_groupName,' vs ',NL_groupName), 
 angle.x = 45, remove.isolate = T,color.heatmap='Spectral')
x11()
gg1 + gg2
dev.copy2pdf(file=paste0(source.cellType,'.sentLRs.',LS_groupName,'.vs.',NL_groupName,'.bubblePlot.pdf'))
#------------------------
NL_groupName <- c('CTRL')
LS_groupName <- c('CAD_without_T2D')

NL_groupName <- c('CTRL')
LS_groupName <- c('CAD_with_T2D')

NL_groupName <- c('CAD_without_T2D')
LS_groupName <- c('CAD_with_T2D')


gg1 <- netVisual_bubble(cellchat, sources.use = levels(cellchat@idents$joint), targets.use = levels(cellchat@idents$joint),  
 comparison = c(dataSetVec[NL_groupName],dataSetVec[LS_groupName]), max.dataset = dataSetVec[LS_groupName],
  title.name = paste0("Increased signaling in ",LS_groupName,' vs ',NL_groupName), 
 angle.x = 45, remove.isolate = T,color.heatmap='Spectral')
gg2 <- netVisual_bubble(cellchat, sources.use = source.cellType, targets.use = target.cellType,
comparison = c(dataSetVec[NL_groupName],dataSetVec[LS_groupName]), max.dataset = dataSetVec[NL_groupName], 
title.name = paste0("Decreased signaling in ",LS_groupName,' vs ',NL_groupName), 
 angle.x = 45, remove.isolate = T,color.heatmap='Spectral')
write.table(gg1$data,file=paste0('All2all.',LS_groupName,'.vs.',NL_groupName,'.increasedSignaling.tsv'),row.names=F,col.names=T,quote=F,sep='\t')
write.table(gg2$data,file=paste0('All2all.',LS_groupName,'.vs.',NL_groupName,'.decreasedSignaling.tsv'),row.names=F,col.names=T,quote=F,sep='\t')


#identify the upgulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs and visulize using bubble plots
#The increased signaling means these signaling have higher communication probability (strength) in one dataset compared to the other dataset.

#------Method2
# identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential gene expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset

pos.dataset = "CAD_without_T2D"
ctrl.dataset = "CTRL_old"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat.CTRL_old.CAD_without_T2D <- identifyOverExpressedGenes(cellchat.CTRL_old.CAD_without_T2D, group.dataset = "datasets", pos.dataset = pos.dataset,
features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# DEG info
cellchat.CTRL_old.CAD_without_T2D@var.features[[paste0(features.name,".info")]]

#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net.CTRL_old.CAD_without_T2D <- netMappingDEG(cellchat.CTRL_old.CAD_without_T2D, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up.CTRL_old.CAD_without_T2D <- subsetCommunication(cellchat.CTRL_old.CAD_without_T2D, net = net.CTRL_old.CAD_without_T2D, 
datasets = pos.dataset,ligand.logFC = 0.1, receptor.logFC = -0.1)
net.up.CTRL_old.CAD_without_T2D <- subset(net.up.CTRL_old.CAD_without_T2D,receptor.logFC!='NA')
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down.CTRL_old.CAD_without_T2D <- subsetCommunication(cellchat.CTRL_old.CAD_without_T2D, net = net.CTRL_old.CAD_without_T2D, 
datasets = ctrl.dataset,ligand.logFC = -0.1, receptor.logFC = -0.1)
net.down.CTRL_old.CAD_without_T2D <- subset(net.down.CTRL_old.CAD_without_T2D,receptor.logFC!='NA')
# Since the signaling genes in the net.up and net.down might be complex with multi-subunits, we can do further deconvolution to obtain the individual signaling genes.
gene.up.CTRL_old.CAD_without_T2D <- extractGeneSubsetFromPair(net.up.CTRL_old.CAD_without_T2D, cellchat.CTRL_old.CAD_without_T2D)
gene.down.CTRL_old.CAD_without_T2D <- extractGeneSubsetFromPair(net.down.CTRL_old.CAD_without_T2D, cellchat.CTRL_old.CAD_without_T2D)
write.table(net.up.CTRL_old.CAD_without_T2D,file='CAD_without_T2D.vs.CTRL_old_net.up.tsv',row.names=F,col.names=T,quote=F,sep='\t')
write.table(net.down.CTRL_old.CAD_without_T2D,file='CAD_without_T2D.vs.CTRL_old_net.down.tsv',row.names=F,col.names=T,quote=F,sep='\t')
#Visualize the enriched ligands, signaling,or ligand-receptor pairs in one condition compared to another condition using wordcloud
# visualize the enriched ligands in the first condition
x11()
computeEnrichmentScore(net.down.CTRL_old.CAD_without_T2D, species = 'human')
dev.copy2pdf(file="CAD_without_T2D.vs.CTRL_old.net.down.wordcloud.pdf")
# visualize the enriched ligands in the second condition
x11()
computeEnrichmentScore(net.up.CTRL_old.CAD_without_T2D, species = 'human')
dev.copy2pdf(file="CAD_without_T2D.vs.CTRL_old.net.up.wordcloud.pdf")



# bubble plot
pairLR.use.up = subset(net.up.CTRL_old.CAD_without_T2D,source %in% c('Adipocyte') & target %in% c('Macrophage','ASPC','Adipocyte'))[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
sources.use = c('Adipocyte'), targets.use = c('Macrophage','ASPC','Adipocyte'),
comparison = c(2, 3),  angle.x = 45, remove.isolate = F,
title.name = paste0("Up-regulated signaling in ", pos.dataset))
gg1$data$interaction_name_2 <- gg1$data$interaction_name
#> Comparing communications on a merged object
pairLR.use.down = subset(net.down.CTRL_old.CAD_without_T2D,source %in% c('Adipocyte') & target %in% c('Macrophage','ASPC','Adipocyte'))[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c('Adipocyte'), targets.use =  c('Macrophage','ASPC','Adipocyte'),
comparison = c(2, 3),  angle.x = 45, remove.isolate = T,
title.name = paste0("Down-regulated signaling in ", pos.dataset))
gg2$data$interaction_name_2 <- gg2$data$interaction_name
#> Comparing communications on a merged object
x11()
gg1 + gg2
dev.copy2pdf(file=paste0(source.cellType,'.sentLRs.',LS_groupName,'.vs.',NL_groupName,'.bubblePlot.pdf'))



#Part IV: Visually compare cell-cell communication using Chord diagram
#Edge color/weight, node color/size/shape: In all visualization plots, edge colors are consistent with the sources as sender, 
#and edge weights are proportional to the interaction strength. Thicker edge line indicates a stronger signal.
# In the Hierarchy plot and Circle plot, circle sizes are proportional to the number of cells in each cell group.
# In the hierarchy plot, solid and open circles represent source and target, respectively.
# In the Chord diagram, the inner thinner bar colors represent the targets that receive signal from the corresponding outer bar. 
# The inner bar size is proportional to the signal strength received by the targets. 
# Such inner bar is helpful for interpreting the complex chord diagram. 
# compare all the interactions sending from fibroblast to inflamatory immune cells

> levels(cellchat@idents$CTRL)
 [1] "Macrophage"  "Adipocyte"   "ASPC"        "TC"          "BC"
 [6] "VEC"         "LEC"         "Neural"      "Mesothelial" "Pericyte"
[11] "NK"          "DC_Neu"      "Plasma"      "SMC"         "Mast"


#  netVisual_chord_gene is used for visualizing the cell-cell communication mediated by mutiple ligand-receptors or signaling pathways 
#(where each sector in the chord diagram is a ligand, receptor or signaling pathway.)
source_cellType = 'Adipocyte'
source_cellType = 'ASPC'
source_cellType = 'Macrophage'

x11()
netVisual_chord_gene(object.list[[1]], sources.use = source_cellType,
  slot.name = "netP",
  legend.pos.x = 10,
  title.name = paste0("Signaling sent from ",source_cellType ,"\n" ,names(object.list)[1]))
dev.copy2pdf(file=paste0(source_cellType, ".sent_signal_pathway",'.CTRL.netVisual_chord_gene.pdf'))


x11()
netVisual_chord_gene(object.list[[2]], sources.use = source_cellType,
  slot.name = "netP",
  legend.pos.x = 10,
  title.name = paste0("Signaling sent from ",source_cellType ,"\n" ,names(object.list)[2]))
dev.copy2pdf(file=paste0(source_cellType, ".sent_signal_pathway",'.CAD_without_T2D.netVisual_chord_gene.pdf'))

x11()
netVisual_chord_gene(object.list[[3]], sources.use = source_cellType,
  slot.name = "netP",
  legend.pos.x = 10,
  title.name = paste0("Signaling sent from ",source_cellType ,"\n" ,names(object.list)[3]))
dev.copy2pdf(file=paste0(source_cellType, ".sent_signal_pathway",'.CAD_with_T2D.netVisual_chord_gene.pdf'))

#  netVisual_chord_gene plot received_signal_pathway
target_cellType = 'Adipocyte'
target_cellType = 'ASPC'
target_cellType = 'Macrophage'

x11()
netVisual_chord_gene(object.list[[1]], targets.use = target_cellType, sources.use = c('Adipocyte','ASPC','Macrophage','TC'),
  slot.name = "netP",
  legend.pos.x = 10,
  title.name = paste0("Signaling received by ",target_cellType ,"\n" ,names(object.list)[1]))
dev.copy2pdf(file=paste0(target_cellType, ".received_signal_pathway",'.CTRL.netVisual_chord_gene.pdf'))

x11()
netVisual_chord_gene(object.list[[2]], targets.use = source_cellType,sources.use = c('Adipocyte','ASPC','Macrophage','TC'),
  slot.name = "netP",
  legend.pos.x = 10,
  title.name = paste0("Signaling received by ",target_cellType ,"\n" ,names(object.list)[2]))
dev.copy2pdf(file=paste0(target_cellType, ".received_signal_pathway",'.CAD_without_T2D.netVisual_chord_gene.pdf'))

x11()
netVisual_chord_gene(object.list[[3]], targets.use = target_cellType,sources.use = c('Adipocyte','ASPC','Macrophage','TC'),
  slot.name = "netP",
  legend.pos.x = 10,
  title.name = paste0("Signaling received by ",target_cellType ,"\n" ,names(object.list)[3]))
dev.copy2pdf(file=paste0(target_cellType, ".received_signal_pathway",'.CAD_with_T2D.netVisual_chord_gene.pdf'))


# Chord diagram
pathways.show <- "SEMA3"
pathways.show <- "VISFATIN"
pathways.show <- "ADIPONECTIN"
pathways.show <- "LEP"
pathways.show <- "IL1"
pathways.show <- "CXCL"
pathways.show <- "CCL"
pathways.show <- "VEGF"
pathways.show <- "TGFb"
pathways.show <- "ANNEXIN"
pathways.show <- "LIFR"
pathways.show <- "EGF"
pathways.show <- "PDGF"
gg1 <- netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "chord", 
signaling.name = paste(pathways.show, names(object.list)[1]))
gg2 <- netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "chord", 
signaling.name = paste(pathways.show, names(object.list)[2]))
gg3 <- netVisual_aggregate(object.list[[3]], signaling = pathways.show, layout = "chord", 
signaling.name = paste(pathways.show, names(object.list)[3]))

x11()
gg1 
dev.copy2pdf(file=paste0(pathways.show,'.CTRL.pathway_chord_diagram.pdf'))
x11()
gg2
dev.copy2pdf(file=paste0(pathways.show,'.CAD_without_T2D.pathway_chord_diagram.pdf'))
x11()
gg3
dev.copy2pdf(file=paste0(pathways.show,'.CAD_with_T2D.pathway_chord_diagram.pdf'))

# Compute the network centrality scores,Visualize the computed centrality scores using heatmap
# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
pathways.show <- "SEMA3"
pathways.show <- "VISFATIN"
pathways.show <- "IL1"
pathways.show <- "CXCL"
pathways.show <- "ANNEXIN"
cellchat.CAD_without_T2D <- netAnalysis_computeCentrality(cellchat.CAD_without_T2D, slot.name = "netP") 
x11()
netAnalysis_signalingRole_network(cellchat.CAD_without_T2D, signaling = pathways.show, width = 15, height = 5, font.size = 10)
dev.copy2pdf(file=paste0(pathways.show,"_CAD_without_T2D_centralityScoresHeatmap.pdf"))

cellchat.CAD_with_T2D <- netAnalysis_computeCentrality(cellchat.CAD_with_T2D, slot.name = "netP") 
x11()
netAnalysis_signalingRole_network(cellchat.CAD_with_T2D, signaling = pathways.show, width = 15, height = 5, font.size = 10)
dev.copy2pdf(file=paste0(pathways.show,"_CAD_with_T2D_centralityScoresHeatmap.pdf"))

cellchat.CTRL <- netAnalysis_computeCentrality(cellchat.CTRL, slot.name = "netP") 
x11()
netAnalysis_signalingRole_network(cellchat.CTRL, signaling = pathways.show, width = 15, height = 5, font.size = 10)
dev.copy2pdf(file=paste0(pathways.show,"_CTRL_centralityScoresHeatmap.pdf"))



#Compute the contribution of each ligand-receptor pair to the overall signaling pathway
pathways.show <- "VISFATIN"
pathways.show <- "SEMA3"
pathways.show <- "ADIPONECTIN"
pathways.show <- "LEP"
pathways.show <- "IL1"
pathways.show <- "CXCL"
pathways.show <- "VEGF"
pathways.show <- "TGFb"
pathways.show <- "ANNEXIN"
pathways.show <- "LIFR"
pathways.show <- "EGF"
pathways.show <- "PDGF"
gg1 <- netAnalysis_contribution(cellchat.CTRL, signaling = pathways.show,title = "Contribution of each L-R pair in CTRL")
gg2 <- netAnalysis_contribution(cellchat.CAD_without_T2D, signaling = pathways.show, title = "Contribution of each L-R pair in CAD_without_T2D")
gg3 <- netAnalysis_contribution(cellchat.CAD_with_T2D, signaling = pathways.show, title = "Contribution of each L-R pair in CAD_with_T2D")
x11()
gg1+gg2+gg3
dev.copy2pdf(file=paste0(pathways.show,"_LRcontribution.pdf"))

pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
pairLR
   interaction_name
1        NAMPT_INSR
2 NAMPT_ITGA5_ITGB1

# LR chord plot

pathways.show = "VISFATIN"
LR.show = 'NAMPT_INSR'

pathways.show = "SEMA3"
LR.show ='SEMA3B_NRP1_PLXNA1'
LR.show ='SEMA3C_NRP1_PLXNA4'
LR.show ='SEMA3A_NRP1_PLXNA4'
LR.show ='SEMA3B_NRP1_PLXNA2'
LR.show ='SEMA3B_NRP1_PLXNA4'

pathways.show <- c("CXCL")
LR.show = 'CXCL12_CXCR4'

pathways.show <- "ANNEXIN"
LR.show = 'ANXA1_FPR1'


pathways.show <- "ANNEXIN"
LR.show = 'ANXA1_FPR1'

pathways.show <- "EGF"
LR.show = 'HBEGF_EGFR'

x11()
netVisual_individual(cellchat.CTRL, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.copy2pdf(file=paste0(pathways.show,".",LR.show,"_chordPlot_CTRL.pdf"))
x11()
netVisual_individual(cellchat.CAD_without_T2D, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.copy2pdf(file=paste0(pathways.show,".",LR.show,"_chordPlot_CAD_without_T2D.pdf"))
x11()
netVisual_individual(cellchat.CAD_with_T2D, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
dev.copy2pdf(file=paste0(pathways.show,".",LR.show,"_chordPlot_CAD_with_T2D.pdf"))


#Part V: Compare the signaling gene expression distribution between different datasets
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
pdf('TGFb.pdf')
plotGeneExpression(cellchat, signaling = "TGFb", split.by = "datasets", colors.ggplot = T)
dev.off()

pdf('PARs.pdf')
plotGeneExpression(cellchat, signaling = "PARs", split.by = "datasets", colors.ggplot = T)
dev.off()

pdf('ANGPTL.pdf')
plotGeneExpression(cellchat, signaling = "ANGPTL", split.by = "datasets", colors.ggplot = T)
dev.off()

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
pdf('ncWNT.pdf')
plotGeneExpression(cellchat, signaling = "ncWNT", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('NOTCH.pdf')
plotGeneExpression(cellchat, signaling = "NOTCH", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('FGF.pdf')
plotGeneExpression(cellchat, signaling = "FGF", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('SEMA3.pdf')
plotGeneExpression(cellchat, signaling = "SEMA3", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('MK.pdf')
plotGeneExpression(cellchat, signaling = "MK", split.by = "datasets", colors.ggplot = T)
dev.off()
pdf('LAMININ.pdf')
plotGeneExpression(cellchat, signaling = "LAMININ", split.by = "datasets", colors.ggplot = T)
dev.off()

pdf('ncWNT.circlePlot.pdf')
pathways.show <- c("ncWNT")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pdf('NOTCH.circlePlot.pdf')
pathways.show <- c("NOTCH")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pdf('FGF.circlePlot.pdf')
pathways.show <- c("FGF")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pdf('SEMA3.circlePlot.pdf')
pathways.show <- c("SEMA3")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


pathways.show <- c("ncWNT")
pairLR.enriched.LS <- extractEnrichedLR(cellchat.LS, signaling = pathways.show, geneLR.return = FALSE)
pairLR.enriched.NL <- extractEnrichedLR(cellchat.NL, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- 'WNT5A_FZD4'
pdf('WNT5A_FZD4.pdf')
for (i in 1:length(object.list)) {
   netVisual_individual(object.list[[i]], signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
}
dev.off()


