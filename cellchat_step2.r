library(CellChat)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)
options(stringsAsFactors = FALSE)
load('cellchat.CAD_with_T2D.ori.Rdata')
load('cellchat.CAD_without_T2D.ori.Rdata')
load('cellchat.CTRL.ori.Rdata')

## Add cell information into meta slot of the object
##If cell mata information is not added when creating CellChat object, USERS can also add it later using addMeta, and set the default cell identities using setIdent.
#cellchat <- addMeta(cellchat, meta = meta)
cellchat.CAD_with_T2D <- setIdent(cellchat.CAD_with_T2D, ident.use = "anno") # set "anno" as default cell identity
cellchat.CAD_without_T2D <- setIdent(cellchat.CAD_without_T2D, ident.use = "anno") # set "anno" as default cell identity
cellchat.CTRL <- setIdent(cellchat.CTRL, ident.use = "anno") # set "anno" as default cell identity


levels(cellchat.CAD_with_T2D@idents) # show factor levels of the cell labels
levels(cellchat.CAD_without_T2D@idents)
levels(cellchat.CTRL@idents)

groupSize.CAD_with_T2D <- as.numeric(table(cellchat.CAD_with_T2D@idents)) # number of cells in each cell group
groupSize.CAD_without_T2D <- as.numeric(table(cellchat.CAD_without_T2D@idents)) # number of cells in each cell group
groupSize.CTRL <- as.numeric(table(cellchat.CTRL@idents)) 

> levels(cellchat.CAD_without_T2D@idents)
 [1] "Macrophage"  "Adipocyte"   "ASPC"        "TC"          "BC"
 [6] "VEC"         "LEC"         "Neural"      "Mesothelial" "Pericyte"
[11] "NK"          "DC_Neu"      "Plasma"      "SMC"         "Mast"

> groupSize.CAD_with_T2D
 [1] 6577 6364 5634 1892  601  863  523 1214  647  596  276  309  163  147   98
> groupSize.CAD_without_T2D
 [1] 7095 4795 5517 1692  308  817  529  787  437  474  255  283  117   93  139
> groupSize.CTRL
 [1] 3694 4734 4686 3098 2295 1382 1432  197 1016  590  313  134  372  118   83

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
#showDatabaseCategory(CellChatDB)

# Show the structure of the database
> dplyr::glimpse(CellChatDB$interaction)
Rows: 1,939
Columns: 11
$ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF…
$ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T…
$ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB…
$ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF…
$ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T…
$ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago…
$ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""…
$ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep…
$ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435…
$ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre…
$ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2…

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.secreted <- subsetDB(CellChatDB, search = "Secreted Signaling",key = "annotation" ) # use Secreted Signaling
> dim(CellChatDB.secreted$interaction)
[1] 1199   11

library(plyr)
#CellChatDB$interaction$anno2 <- mapvalues(CellChatDB$interaction$annotation,from=c('Secreted Signaling','Cell-Cell Contact'),to=c("Considered","Considered"))
#CellChatDB.use <- subsetDB(CellChatDB, search = "Considered",key = "anno2" )
#dplyr::glimpse(CellChatDB.use$interaction)
#> dplyr::glimpse(CellChatDB.use$interaction)
#Rows: 1,518
#Columns: 12

# use all CellChatDB for cell-cell communication analysis
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB (all pathways)
CellChatDB.use <- CellChatDB.secreted # simply use Secreted Signaling CellChatDB
# set the used database in the object
cellchat.CAD_with_T2D@DB <- CellChatDB.use
cellchat.CAD_without_T2D@DB <- CellChatDB.use
cellchat.CTRL@DB <- CellChatDB.use


# Preprocessing the expression data for cell-cell communication analysis
#To infer the cell state-specific communications, we identify over-expressed ligands or receptors in one cell group and
#then identify over-expressed ligand-receptor interactions if either ligand or receptor is over-expressed.
#Project gene expression data onto protein-protein interaction (PPI) network
# subset the expression data of signaling genes for saving computation cost
#An updated CellChat object by assigning a subset of the data into the slot `data.signaling`
cellchat.CAD_with_T2D <- subsetData(cellchat.CAD_with_T2D)
cellchat.CAD_without_T2D <- subsetData(cellchat.CAD_without_T2D) # subset the expression data of signaling genes for saving computation cost
cellchat.CTRL <- subsetData(cellchat.CTRL)

future::plan("multiprocess", workers = 50) # do parallel
cellchat.CAD_with_T2D <- identifyOverExpressedGenes(cellchat.CAD_with_T2D)
cellchat.CAD_with_T2D <- identifyOverExpressedInteractions(cellchat.CAD_with_T2D)
cellchat.CAD_with_T2D <- projectData(cellchat.CAD_with_T2D, PPI.human) #PPI.human for human

cellchat.CAD_without_T2D <- identifyOverExpressedGenes(cellchat.CAD_without_T2D)
cellchat.CAD_without_T2D <- identifyOverExpressedInteractions(cellchat.CAD_without_T2D)
cellchat.CAD_without_T2D <- projectData(cellchat.CAD_without_T2D, PPI.human) #PPI.human for human

cellchat.CTRL <- identifyOverExpressedGenes(cellchat.CTRL)
cellchat.CTRL <- identifyOverExpressedInteractions(cellchat.CTRL)
cellchat.CTRL <- projectData(cellchat.CTRL, PPI.human) #PPI.human for human


## Part II: Inference of cell-cell communication network

#CellChat infers the biologically significant cell-cell communication by assigning each interaction with a probability value
#and peforming a permutation test. CellChat models the probability of cell-cell communication by integrating gene expression
#with prior known knowledge of the interactions between signaling ligands, receptors and their cofactors using the law of mass action.
# Compute the communication probability and infer cellular communication network

#Set raw.use= FALSE to use the projected data when analyzing single-cell
#          data with shallow sequencing depth because the projected data
#          could help to reduce the dropout effects of signaling genes,
#          in particular for possible zero expression of subunits of
#          ligands/receptors.
#Set population.size =
#          TRUE if analyzing unsorted single-cell transcriptomes, with
#          the reason that abundant cell populations tend to send
#         collectively stronger signals than the rare cell populations.

#CellChat uses a statistically robust mean method called ‘trimean’, which produces fewer interactions than other methods. 
#However, we find that CellChat performs well at### predicting stronger interactions,
# which is very helpful for narrowing down on interactions for further experimental validations. 
#In computeCommunProb, we provide an option for using other methods, such as 5% and 10% truncated mean, 
#to calculating the average gene expression. Of note, ‘trimean’ approximates 25% truncated mean,
# implying that the average gene expression is zero if the percent of expressed cells in one group is less than 25%.

#defalt setting
#cellchat.NL <- computeCommunProb(cellchat.NL, raw.use = TRUE,population.size=FALSE) #trimean
#cellchat.LS <- computeCommunProb(cellchat.LS, raw.use = TRUE,population.size=FALSE) #trimean
#cellchat.NL.trimean <- computeCommunProb(cellchat.NL, raw.use = FALSE,population.size=FALSE) #trimean
#cellchat.LS.trimean <- computeCommunProb(cellchat.LS, raw.use = FALSE,population.size=FALSE) #trimean
#save(cellchat.NL.trimean,file='cellchat.NL.trimean.Rdata')
#save(cellchat.LS.trimean,file='cellchat.LS.trimean.Rdata')

cellchat.CAD_without_T2D <- computeCommunProb(cellchat.CAD_without_T2D, raw.use = TRUE,population.size=TRUE,type = "truncatedMean",trim = 0.1)#10% truncated mean; to get more interactions
cellchat.CAD_with_T2D <- computeCommunProb(cellchat.CAD_with_T2D, raw.use = TRUE,population.size=TRUE,type = "truncatedMean",trim = 0.1)#10% truncated mean
cellchat.CTRL <- computeCommunProb(cellchat.CTRL, raw.use = TRUE,population.size=TRUE,type = "truncatedMean",trim = 0.1)#10% truncated mean


# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat.CAD_without_T2D <- filterCommunication(cellchat.CAD_without_T2D, min.cells = 10)
cellchat.CAD_with_T2D <- filterCommunication(cellchat.CAD_with_T2D, min.cells = 10)
cellchat.CTRL <- filterCommunication(cellchat.CTRL, min.cells = 10)


#Extract the inferred cellular communication network as a data frame
df.net.CAD_without_T2D <- subsetCommunication(cellchat.CAD_without_T2D)
df.net.CAD_with_T2D <- subsetCommunication(cellchat.CAD_with_T2D)
df.net.CTRL <- subsetCommunication(cellchat.CTRL)

write.table(df.net.CAD_without_T2D,file='df.net.CAD_without_T2D.tsv',row.names=F,col.name=T,sep='\t',quote=F)
write.table(df.net.CAD_with_T2D,file='df.net.CAD_with_T2D.tsv',row.names=F,col.name=T,sep='\t',quote=F)
write.table(df.net.CTRL,file='df.net.CTRL.tsv',row.names=F,col.name=T,sep='\t',quote=F)


#Infer the cell-cell communication at a signaling pathway level
#CellChat computes the communication probability on signaling pathway level by summarizing the communication probabilities of all ligands-receptors interactions associated with each signaling pathway.
#The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
cellchat.CAD_without_T2D <- computeCommunProbPathway(cellchat.CAD_without_T2D)
cellchat.CAD_with_T2D <- computeCommunProbPathway(cellchat.CAD_with_T2D)
cellchat.CTRL <- computeCommunProbPathway(cellchat.CTRL)



#Calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability
cellchat.CAD_without_T2D <- aggregateNet(cellchat.CAD_without_T2D)
cellchat.CAD_with_T2D <- aggregateNet(cellchat.CAD_with_T2D)
cellchat.CTRL <- aggregateNet(cellchat.CTRL)


save(cellchat.CAD_without_T2D,file='cellchat.CAD_without_T2D.RData')
save(cellchat.CAD_with_T2D,file='cellchat.CAD_with_T2D.RData')
save(cellchat.CTRL,file='cellchat.CTRL.RData')




#visualize the aggregated cell-cell communication network. For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
x11()
netVisual_circle(cellchat.CAD_without_T2D@net$count, vertex.weight = groupSize.CAD_without_T2D, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.copy2pdf(file='cellchat.CAD_without_T2D.network.NumInteraction.pdf')
x11()
netVisual_circle(cellchat.CAD_with_T2D@net$count, vertex.weight = groupSize.CAD_with_T2D, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.copy2pdf(file='cellchat.CAD_with_T2D.network.NumInteraction.pdf')

x11()
netVisual_circle(cellchat.CTRL@net$count, vertex.weight = groupSize.CTRL, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.copy2pdf(file='cellchat.CTRL.network.NumInteraction.pdf')



pdf(file='cellchat.CTRL.network.NumInteraction.eachGroup.pdf',width=10,height=20)
mat.CTRL <- cellchat.CTRL@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat.CTRL)) {
  mat2 <- matrix(0, nrow = nrow(mat.CTRL), ncol = ncol(mat.CTRL), dimnames = dimnames(mat.CTRL))
  mat2[i, ] <- mat.CTRL[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize.CTRL, weight.scale = T, edge.weight.max = max(mat.CTRL), title.name = rownames(mat.CTRL)[i])
}
dev.off()


pdf(file='cellchat.CAD_with_T2D.network.NumInteraction.eachGroup.pdf',width=10,height=20)
mat.CAD_with_T2D <- cellchat.CAD_with_T2D@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat.CAD_with_T2D)) {
  mat2 <- matrix(0, nrow = nrow(mat.CAD_with_T2D), ncol = ncol(mat.CAD_with_T2D), dimnames = dimnames(mat.CAD_with_T2D))
  mat2[i, ] <- mat.CAD_with_T2D[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize.CAD_with_T2D, weight.scale = T, edge.weight.max = max(mat.CAD_with_T2D), title.name = rownames(mat.CAD_with_T2D)[i])
}
dev.off()

pdf(file='cellchat.CAD_without_T2D.network.NumInteraction.eachGroup.pdf',width=10,height=20)
mat.CAD_without_T2D <- cellchat.CAD_without_T2D@net$weight
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat.CAD_without_T2D)) {
  mat2 <- matrix(0, nrow = nrow(mat.CAD_without_T2D), ncol = ncol(mat.CAD_without_T2D), dimnames = dimnames(mat.CAD_without_T2D))
  mat2[i, ] <- mat.CAD_without_T2D[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize.CAD_without_T2D, weight.scale = T, edge.weight.max = max(mat.CAD_without_T2D), title.name = rownames(mat.CAD_without_T2D)[i])
}
dev.off()


#Identify global communication patterns to explore how multiple cell types and signaling pathways coordinate together
#CellChat employs a pattern recognition method to identify the global communication patterns.
#river plots
library(NMF)
library(ggalluvial)
#Here we run selectK to infer the number of patterns.

############### CAD_with_T2D
#selectK(cellchat.CAD_with_T2D, pattern = "outgoing")
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 7.
nPatterns = 15
cellchat.CAD_with_T2D <- identifyCommunicationPatterns(cellchat.CAD_with_T2D, pattern = "outgoing", k = nPatterns)
# river plot
x11()
netAnalysis_river(cellchat.CAD_with_T2D, pattern = "outgoing")
dev.copy2pdf(file="CAD_with_T2D.riverPlot.outgoingPattern.pdf")
# dot plot
x11()
netAnalysis_dot(cellchat.CAD_with_T2D, pattern = "outgoing")
dev.copy2pdf(file="CAD_with_T2D.dotPlot.outgoingPattern.pdf")

#selectK(cellchat.CAD_with_T2D, pattern = "incoming")
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 8.
nPatterns = 15
cellchat.CAD_with_T2D <- identifyCommunicationPatterns(cellchat.CAD_with_T2D, pattern = "incoming", k = nPatterns)
# river plot
x11()
netAnalysis_river(cellchat.CAD_with_T2D, pattern = "incoming")
dev.copy2pdf(file="CAD_with_T2D.riverPlot.incomingPattern.pdf")
# dot plot
x11()
netAnalysis_dot(cellchat.CAD_with_T2D, pattern = "incoming")
dev.copy2pdf(file="CAD_with_T2D.dotPlot.incomingPattern.pdf")

############### CAD_without_T2D
#selectK(cellchat.CAD_without_T2D, pattern = "outgoing")
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 7.
nPatterns = 15
cellchat.CAD_without_T2D <- identifyCommunicationPatterns(cellchat.CAD_without_T2D, pattern = "outgoing", k = nPatterns)
# river plot
x11()
netAnalysis_river(cellchat.CAD_without_T2D, pattern = "outgoing")
dev.copy2pdf(file="CAD_without_T2D.riverPlot.outgoingPattern.pdf")
# dot plot
x11()
netAnalysis_dot(cellchat.CAD_without_T2D, pattern = "outgoing")
dev.copy2pdf(file="CAD_without_T2D.dotPlot.outgoingPattern.pdf")

#selectK(cellchat.CAD_without_T2D, pattern = "incoming")
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 8.
nPatterns = 15
cellchat.CAD_without_T2D <- identifyCommunicationPatterns(cellchat.CAD_without_T2D, pattern = "incoming", k = nPatterns)
# river plot
x11()
netAnalysis_river(cellchat.CAD_without_T2D, pattern = "incoming")
dev.copy2pdf(file="CAD_without_T2D.riverPlot.incomingPattern.pdf")
# dot plot
x11()
netAnalysis_dot(cellchat.CAD_without_T2D, pattern = "incoming")
dev.copy2pdf(file="CAD_without_T2D.dotPlot.incomingPattern.pdf")


############### CTRL
#selectK(cellchat.CTRL, pattern = "outgoing")
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 7.
nPatterns = 15
cellchat.CTRL <- identifyCommunicationPatterns(cellchat.CTRL, pattern = "outgoing", k = nPatterns)
# river plot
x11()
netAnalysis_river(cellchat.CTRL, pattern = "outgoing")
dev.copy2pdf(file="CTRL.riverPlot.outgoingPattern.pdf")
# dot plot
x11()
netAnalysis_dot(cellchat.CTRL, pattern = "outgoing")
dev.copy2pdf(file="CTRL.dotPlot.outgoingPattern.pdf")

#selectK(cellchat.CTRL, pattern = "incoming")
#Both Cophenetic and Silhouette values begin to drop suddenly when the number of outgoing patterns is 8.
nPatterns = 15
cellchat.CTRL <- identifyCommunicationPatterns(cellchat.CTRL, pattern = "incoming", k = nPatterns)
# river plot
x11()
netAnalysis_river(cellchat.CTRL, pattern = "incoming")
dev.copy2pdf(file="CTRL.riverPlot.incomingPattern.pdf")
# dot plot
x11()
netAnalysis_dot(cellchat.CTRL, pattern = "incoming")
dev.copy2pdf(file="CTRL.dotPlot.incomingPattern.pdf")