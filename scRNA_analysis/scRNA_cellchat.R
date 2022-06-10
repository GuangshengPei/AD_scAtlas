library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)

# Create a CellChat object
# C1:PFC_Health C2:EC_Health C3:SFG_Health, D1:PFC_Patient D2:EC_Patient D3:cellchat_D3
cellchat_C1 <- createCellChat(object = data.input_C1, meta = identity_C1, group.by = "labels")
cellchat_C2 <- createCellChat(object = data.input_C2, meta = identity_C2, group.by = "labels")
cellchat_C3 <- createCellChat(object = data.input_C3, meta = identity_C3, group.by = "labels")
cellchat_D1 <- createCellChat(object = data.input_D1, meta = identity_D1, group.by = "labels")
cellchat_D2 <- createCellChat(object = data.input_D2, meta = identity_D2, group.by = "labels")
cellchat_D3 <- createCellChat(object = data.input_D3, meta = identity_D3, group.by = "labels")

#C1
cellchat_C1 <- addMeta(cellchat_C1, meta = identity_C1)
cellchat_C1 <- setIdent(cellchat_C1, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_C1@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_C1@idents)) # number of cells in each cell group

cellchat_C1DB <- CellChatDB.human # 
showDatabaseCategory(cellchat_C1DB)
colnames(cellchat_C1DB$interaction)

cellchat_C1DB.use <- cellchat_C1DB
cellchat_C1DB.use <- subsetDB(cellchat_C1DB, search = "Secreted Signaling")
cellchat_C1@DB <- cellchat_C1DB.use # set the used database in the object

unique(cellchat_C1DB$interaction$annotation)
cellchat_C1 <- subsetData(cellchat_C1) # subset the expression data of signaling genes for saving computation cost
#Preprocessing the expression data for cell-cell communication analysis
future::plan("multiprocess", workers = 10) # do parallel
cellchat_C1 <- identifyOverExpressedGenes(cellchat_C1)
cellchat_C1 <- identifyOverExpressedInteractions(cellchat_C1)
cellchat_C1 <- projectData(cellchat_C1, PPI.human)

cellchat_C1 <- computeCommunProb(cellchat_C1, raw.use = TRUE)
cellchat_C1 <- filterCommunication(cellchat_C1, min.cells = 10)
cellchat_C1 <- computeCommunProbPathway(cellchat_C1)
cellchat_C1 <- aggregateNet(cellchat_C1)

#D1
cellchat_D1 <- addMeta(cellchat_D1, meta = identity_D1)
cellchat_D1 <- setIdent(cellchat_D1, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_D1@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_D1@idents)) # number of cells in each cell group

cellchat_D1DB <- CellChatDB.human # 
showDatabaseCategory(cellchat_D1DB)
colnames(cellchat_D1DB$interaction)

cellchat_D1DB.use <- cellchat_D1DB
cellchat_D1DB.use <- subsetDB(cellchat_D1DB, search = "Secreted Signaling")
cellchat_D1@DB <- cellchat_D1DB.use # set the used database in the object

unique(cellchat_D1DB$interaction$annotation)
cellchat_D1 <- subsetData(cellchat_D1) # subset the expression data of signaling genes for saving computation cost
#Preprocessing the expression data for cell-cell communication analysis
future::plan("multiprocess", workers = 10) # do parallel
cellchat_D1 <- identifyOverExpressedGenes(cellchat_D1)
cellchat_D1 <- identifyOverExpressedInteractions(cellchat_D1)
cellchat_D1 <- projectData(cellchat_D1, PPI.human)

cellchat_D1 <- computeCommunProb(cellchat_D1, raw.use = TRUE)
cellchat_D1 <- filterCommunication(cellchat_D1, min.cells = 10)
cellchat_D1 <- computeCommunProbPathway(cellchat_D1)
cellchat_D1 <- aggregateNet(cellchat_D1)

#C2
cellchat_C2 <- addMeta(cellchat_C2, meta = identity_C2)
cellchat_C2 <- setIdent(cellchat_C2, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_C2@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_C2@idents)) # number of cells in each cell group

cellchat_C2DB <- CellChatDB.human # 
showDatabaseCategory(cellchat_C2DB)
colnames(cellchat_C2DB$interaction)

cellchat_C2DB.use <- cellchat_C2DB
cellchat_C2DB.use <- subsetDB(cellchat_C2DB, search = "Secreted Signaling")
cellchat_C2@DB <- cellchat_C2DB.use # set the used database in the object

unique(cellchat_C2DB$interaction$annotation)
cellchat_C2 <- subsetData(cellchat_C2) # subset the expression data of signaling genes for saving computation cost
#Preprocessing the expression data for cell-cell communication analysis
future::plan("multiprocess", workers = 10) # do parallel
cellchat_C2 <- identifyOverExpressedGenes(cellchat_C2)
cellchat_C2 <- identifyOverExpressedInteractions(cellchat_C2)
cellchat_C2 <- projectData(cellchat_C2, PPI.human)

cellchat_C2 <- computeCommunProb(cellchat_C2, raw.use = TRUE)
cellchat_C2 <- filterCommunication(cellchat_C2, min.cells = 10)
cellchat_C2 <- computeCommunProbPathway(cellchat_C2)
cellchat_C2 <- aggregateNet(cellchat_C2)

#D2
cellchat_D2 <- addMeta(cellchat_D2, meta = identity_D2)
cellchat_D2 <- setIdent(cellchat_D2, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_D2@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_D2@idents)) # number of cells in each cell group

cellchat_D2DB <- CellChatDB.human # 
showDatabaseCategory(cellchat_D2DB)
colnames(cellchat_D2DB$interaction)

cellchat_D2DB.use <- cellchat_D2DB
cellchat_D2DB.use <- subsetDB(cellchat_D2DB, search = "Secreted Signaling")
cellchat_D2@DB <- cellchat_D2DB.use # set the used database in the object

unique(cellchat_D2DB$interaction$annotation)
cellchat_D2 <- subsetData(cellchat_D2) # subset the expression data of signaling genes for saving computation cost
#Preprocessing the expression data for cell-cell communication analysis
future::plan("multiprocess", workers = 10) # do parallel
cellchat_D2 <- identifyOverExpressedGenes(cellchat_D2)
cellchat_D2 <- identifyOverExpressedInteractions(cellchat_D2)
cellchat_D2 <- projectData(cellchat_D2, PPI.human)

cellchat_D2 <- computeCommunProb(cellchat_D2, raw.use = TRUE)
cellchat_D2 <- filterCommunication(cellchat_D2, min.cells = 10)
cellchat_D2 <- computeCommunProbPathway(cellchat_D2)
cellchat_D2 <- aggregateNet(cellchat_D2)

#C3
cellchat_C3 <- addMeta(cellchat_C3, meta = identity_C3)
cellchat_C3 <- setIdent(cellchat_C3, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_C3@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_C3@idents)) # number of cells in each cell group

cellchat_C3DB <- CellChatDB.human # 
showDatabaseCategory(cellchat_C3DB)
colnames(cellchat_C3DB$interaction)

cellchat_C3DB.use <- cellchat_C3DB
cellchat_C3DB.use <- subsetDB(cellchat_C3DB, search = "Secreted Signaling")
cellchat_C3@DB <- cellchat_C3DB.use # set the used database in the object

unique(cellchat_C3DB$interaction$annotation)
cellchat_C3 <- subsetData(cellchat_C3) # subset the expression data of signaling genes for saving computation cost
#Preprocessing the expression data for cell-cell communication analysis
future::plan("multiprocess", workers = 10) # do parallel
cellchat_C3 <- identifyOverExpressedGenes(cellchat_C3)
cellchat_C3 <- identifyOverExpressedInteractions(cellchat_C3)
cellchat_C3 <- projectData(cellchat_C3, PPI.human)

cellchat_C3 <- computeCommunProb(cellchat_C3, raw.use = TRUE)
cellchat_C3 <- filterCommunication(cellchat_C3, min.cells = 10)
cellchat_C3 <- computeCommunProbPathway(cellchat_C3)
cellchat_C3 <- aggregateNet(cellchat_C3)

#D3
cellchat_D3 <- addMeta(cellchat_D3, meta = identity_D3)
cellchat_D3 <- setIdent(cellchat_D3, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_D3@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat_D3@idents)) # number of cells in each cell group

cellchat_D3DB <- CellChatDB.human # 
showDatabaseCategory(cellchat_D3DB)
colnames(cellchat_D3DB$interaction)

cellchat_D3DB.use <- cellchat_D3DB
cellchat_D3DB.use <- subsetDB(cellchat_D3DB, search = "Secreted Signaling")
cellchat_D3@DB <- cellchat_D3DB.use # set the used database in the object

unique(cellchat_D3DB$interaction$annotation)
cellchat_D3 <- subsetData(cellchat_D3) # subset the expression data of signaling genes for saving computation cost
#Preprocessing the expression data for cell-cell communication analysis
future::plan("multiprocess", workers = 20) # do parallel
cellchat_D3 <- identifyOverExpressedGenes(cellchat_D3)
cellchat_D3 <- identifyOverExpressedInteractions(cellchat_D3)
cellchat_D3 <- projectData(cellchat_D3, PPI.human)

cellchat_D3 <- computeCommunProb(cellchat_D3, raw.use = TRUE)
cellchat_D3 <- filterCommunication(cellchat_D3, min.cells = 10)
cellchat_D3 <- computeCommunProbPathway(cellchat_D3)
cellchat_D3 <- aggregateNet(cellchat_D3)

# Can back to laptop 
cellchat_C1 <- netAnalysis_computeCentrality(cellchat_C1, slot.name = "netP")
cellchat_C2 <- netAnalysis_computeCentrality(cellchat_C2, slot.name = "netP")
cellchat_C3 <- netAnalysis_computeCentrality(cellchat_C3, slot.name = "netP")
cellchat_D1 <- netAnalysis_computeCentrality(cellchat_D1, slot.name = "netP")
cellchat_D2 <- netAnalysis_computeCentrality(cellchat_D2, slot.name = "netP")
cellchat_D3 <- netAnalysis_computeCentrality(cellchat_D3, slot.name = "netP")

object.list <- list(PFC_Health = cellchat_C1, EC_Health = cellchat_C2, SFG_Health = cellchat_C3, PFC_Patient = cellchat_D1, EC_Patient = cellchat_D2, SFG_Patient = cellchat_D3)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,4,5,6))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,1,1,2,2,2), measure = "weight")
gg3 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,1,2,3))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,1,1,2,2,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,1,1,2,2,2), measure = "weight")
gg3 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,1,2,3))
gg4 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3,1,2,3), measure = "weight")

pdf("Overview.pdf", 9.5, 4)
gg1 + gg2
gg3 + gg4
dev.off()

gg1 <- netVisual_heatmap(cellchat, comparison = c(1, 4))
gg2 <- netVisual_heatmap(cellchat, comparison = c(2, 5))
gg3 <- netVisual_heatmap(cellchat, comparison = c(3, 6))
gg1 + gg2 + gg3
#> Do heatmap based on a merged object
gg4 <- netVisual_heatmap(cellchat, comparison = c(1, 4), measure = "weight")
gg5 <- netVisual_heatmap(cellchat, comparison = c(2, 5), measure = "weight")
gg6 <- netVisual_heatmap(cellchat, comparison = c(3, 6), measure = "weight")

gg1 + gg2 + gg3
gg4 + gg5 + gg6
#=====================
#The differential network analysis only works for pairwise datasets. If there are more datasets for comparison, we can directly show the number of interactions or interaction strength between any two cell populations in each dataset.
weight.max <- getMaxWeight(object.list, attribute = c("idents", "count"))

pdf("Circle_network.pdf", 12, 6)
par(mfrow = c(2,3), xpd=TRUE)
par(mar=c(2, 1, 2, 1))
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

pdf("_Circle_network_L.pdf", 12, 6)
par(mfrow = c(2,3), xpd=TRUE)
par(mar=c(2, 1, 2, 1))
for (i in 1:length(object.list)) {
	netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

#===============================================================
#Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
pdf("Multi_group_2D.pdf", 12, 6)
patchwork::wrap_plots(plots = gg)
dev.off()

par(mfrow = c(2,3), xpd=TRUE)
pdf("netVisual_diffInteraction.pdf", 7, 7)
par(mar=c(2, 1, 2, 1))
	netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(1, 4), label.edge= T)
	netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(2, 5), label.edge= T)
	netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(3, 6), label.edge= T)
	netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(1, 4), label.edge= T)
	netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(2, 5), label.edge= T)
	netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(3, 6), label.edge= T)
dev.off()

par(mfrow = c(2,3), xpd=TRUE)
pdf("netVisual_diffInteraction_count_weight.pdf", 7, 7)
par(mfrow = c(2,3), xpd=TRUE)
par(mar=c(2, 1, 2, 1))
	netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(1, 4), label.edge= T)
	netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(2, 5), label.edge= T)
	netVisual_diffInteraction(cellchat, weight.scale = T, comparison = c(3, 6), label.edge= T)
	netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(1, 4), label.edge= T)
	netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(2, 5), label.edge= T)
	netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", comparison = c(3, 6), label.edge= T)
dev.off()

#============================================================================
#Part II: Identify the conserved and context-specific signaling pathways
#=======================================
#Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
pdf("NMF_functional.pdf", 8, 7)
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3, point.shape = c(21, 22, 24, 16, 15, 17))
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3, comparison2 = c(1, 4))
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3, comparison = c(2, 5))
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3, comparison = c(3, 6))
netVisual_embeddingPairwiseZoomIn(cellchat, type = "functional", nCol = 2, point.shape = c(21, 22, 24, 16, 15, 17))
dev.off()

#netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

# Long times
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
pdf("NMF_structural.pdf", 8, 7)
netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3, point.shape = c(21, 22, 24, 16, 15, 17))
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2, point.shape = c(21, 22, 24, 16, 15, 17))
dev.off()

#=============================================================================================
#Compute and visualize the pathway distance in the learned joint manifold
#rankSimilarity(cellchat, type = "functional")
library(RColorBrewer)
mycol <- c(brewer.pal(9,"Set1"),brewer.pal(8, "Set2"),brewer.pal(12,"Set3")[c(-2,-12)],brewer.pal(12,"Paired"))

#Identify and visualize the conserved and context-specific signaling pathways
gg1 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = FALSE, comparison = c(1, 4), cutoff.pvalue = 0.05)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = FALSE, comparison = c(2, 5), cutoff.pvalue = 0.05)
gg3 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = FALSE, comparison = c(3, 6), cutoff.pvalue = 0.05)
gg4 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = FALSE, comparison = c(1,2,3,4,5,6))
gg5 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = FALSE, comparison = c(1,4,2,5,3,6))	

pdf("Part 2_Figure 3_rankNet_1.pdf", 15, 8)
gg1 + gg2 + gg3
gg4
dev.off()

pdf("Part 2_Figure 3_rankNet_1_all.pdf", 6, 8)
gg4
gg5
dev.off()

#===
gg6 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = FALSE, comparison = c(1, 4), cutoff.pvalue = 0.05)
gg7 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = FALSE, comparison = c(2, 5), cutoff.pvalue = 0.05)
gg8 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = FALSE, comparison = c(3, 6), cutoff.pvalue = 0.05)
gg9 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = FALSE, comparison = c(1,2,3,4,5,6), cutoff.pvalue = 0.05)
gg10 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = FALSE, comparison = c(1,4,2,5,3,6), cutoff.pvalue = 0.05)

pdf("Part 2_Figure 4_rankNet_2.pdf", 15, 8)
gg6 + gg7 + gg8
gg9
dev.off()

pdf("Part 2_Figure 4_rankNet_2_all.pdf", 6, 8)
gg9
gg10
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(1, 4), return.data = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(2, 5), return.data = TRUE)
gg3 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(3, 6), return.data = TRUE)

gg_all <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = FALSE, comparison = c(1,4,2,5,3,6), return.data = TRUE)
target_signal <- names(which(tapply(gg_all$signaling.contribution$contribution.scaled, gg_all$signaling.contribution$name, max) > 1e-3))
gg_all_clean <- gg_all$signaling.contribution[as.character(gg_all$signaling.contribution$name) %in% target_signal,]

write.table(gg_all_clean, file = "target_signal_all.txt", quote = F, sep = "\t")

target_signal1 <- names(which(tapply(gg1$signaling.contribution$contribution.scaled, gg1$signaling.contribution$name, max) > 1e-3))
target_signal2 <- names(which(tapply(gg2$signaling.contribution$contribution.scaled, gg2$signaling.contribution$name, max) > 1e-3))
target_signal3 <- names(which(tapply(gg3$signaling.contribution$contribution.scaled, gg3$signaling.contribution$name, max) > 1e-3))

gg1_clean <- gg1$signaling.contribution[as.character(gg1$signaling.contribution$name) %in% target_signal1,]
gg2_clean <- gg2$signaling.contribution[as.character(gg2$signaling.contribution$name) %in% target_signal2,]
gg3_clean <- gg3$signaling.contribution[as.character(gg3$signaling.contribution$name) %in% target_signal3,]

gg1_clean$pvalues_strand <- gg1_clean$pvalues
gg1_clean$pvalues_strand[gg1_clean$contribution.relative.1 < 1] <- -1 * gg1_clean$pvalues_strand[gg1_clean$contribution.relative.1 < 1]

gg2_clean$pvalues_strand <- gg2_clean$pvalues
gg2_clean$pvalues_strand[gg2_clean$contribution.relative.1 < 1] <- -1 * gg2_clean$pvalues_strand[gg2_clean$contribution.relative.1 < 1]

gg3_clean$pvalues_strand <- gg3_clean$pvalues
gg3_clean$pvalues_strand[gg3_clean$contribution.relative.1 < 1] <- -1 * gg3_clean$pvalues_strand[gg3_clean$contribution.relative.1 < 1]

write.table(gg1_clean, file = "target_signal1.txt", quote = F, sep = "\t")
write.table(gg2_clean, file = "target_signal2.txt", quote = F, sep = "\t")
write.table(gg3_clean, file = "target_signal3.txt", quote = F, sep = "\t")

signal_1 <- as.character(unique(gg1_clean[gg1_clean$pvalues < 0.05,]$name))
signal_2 <- as.character(unique(gg2_clean[gg2_clean$pvalues < 0.05,]$name))
signal_3 <- as.character(unique(gg3_clean[gg3_clean$pvalues < 0.05,]$name))

signal_all <- data.frame(table(c(signal_1, signal_2, signal_3)))
signal_all <- signal_all[order(signal_all$Freq, decreasing = T),]

#===========================================
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+3]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6, color.heatmap = "GnBu")

pdf("communication_heatmap_region1.pdf", 6.5, 6)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
draw(ht3 + ht4, ht_gap = unit(1.5, "cm"))
dev.off()

i = 2
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+3]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6, color.heatmap = "GnBu")

pdf("communication_heatmap_region2.pdf", 6.5, 6)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
draw(ht3 + ht4, ht_gap = unit(1.5, "cm"))
dev.off()

i = 3
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+3]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6, color.heatmap = "GnBu")

pdf("communication_heatmap_region3.pdf", 6.5, 6)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
draw(ht3 + ht4, ht_gap = unit(1.5, "cm"))
dev.off()

i = 1
pathway.union <- unique(c(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways, object.list[[3]]@netP$pathways, object.list[[4]]@netP$pathways, object.list[[5]]@netP$pathways, object.list[[6]]@netP$pathways))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6, color.heatmap = "GnBu")

pdf("communication_heatmap_region1b.pdf", 6.5, 6)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
draw(ht3 + ht4, ht_gap = unit(1.5, "cm"))
dev.off()

i = 2
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6, color.heatmap = "GnBu")

pdf("communication_heatmap_region2b.pdf", 6.5, 6)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
draw(ht3 + ht4, ht_gap = unit(1.5, "cm"))
dev.off()

i = 3
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6)
ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+3]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+3], width = 5, height = 6, color.heatmap = "GnBu")

pdf("communication_heatmap_region3b.pdf", 6.5, 6)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"))
draw(ht3 + ht4, ht_gap = unit(1.5, "cm"))
dev.off()

#===========================================
#Identify signaling roles (e.g., dominant senders, receivers) of cell groups as well as the major contributing signaling
pathways.show <- c("PTPRM", "NEGR", "NCAM", "CLDN", "CADM", "NRXN", "NRG")

pdf("network centrality_C1.pdf", 8, 6)
netAnalysis_signalingRole_network(cellchat_C1, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()
pdf("network centrality_C2.pdf", 8, 6)
netAnalysis_signalingRole_network(cellchat_C2, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()
pdf("network centrality_C3.pdf", 8, 6)
netAnalysis_signalingRole_network(cellchat_C3, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()

pdf("centrality_D1.pdf", 8, 6)
netAnalysis_signalingRole_network(cellchat_D1, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()
pdf("network centrality_D2.pdf", 8, 6)
netAnalysis_signalingRole_network(cellchat_D2, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()
pdf("network centrality_D3.pdf", 8, 6)
netAnalysis_signalingRole_network(cellchat_D3, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
dev.off()

pdf("netVisual_bubble_1.pdf", 10, 9)
for (i in 1:6){
gg1 <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:6),  comparison = c(1, 4), max.dataset = 4, title.name = "Increased signaling in ALZ PFC", angle.x = 90, remove.isolate = F, thresh = 0.1)
gg2 <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:6),  comparison = c(1, 4), max.dataset = 1, title.name = "Decreased signaling in ALZ PFC", angle.x = 90, remove.isolate = F, thresh = 0.1)
print (gg1 + gg2)
}
dev.off()

pdf("netVisual_bubble_1b.pdf", 10, 9)
for (i in 1:6){
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = i,  comparison = c(1, 4), max.dataset = 4, title.name = "Increased signaling in ALZ PFC", angle.x = 90, remove.isolate = F, thresh = 0.1)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = i,  comparison = c(1, 4), max.dataset = 1, title.name = "Decreased signaling in ALZ PFC", angle.x = 90, remove.isolate = F, thresh = 0.1)
print (gg1 + gg2)
}
dev.off()

pdf("netVisual_bubble_2.pdf", 10, 9)
for (i in 1:6){
gg1 <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:6),  comparison = c(2, 5), max.dataset = 5, title.name = "Increased signaling in ALZ EC", angle.x = 90, remove.isolate = F, thresh = 0.1)
gg2 <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:6),  comparison = c(2, 5), max.dataset = 2, title.name = "Decreased signaling in ALZ EC", angle.x = 90, remove.isolate = F, thresh = 0.1)
print (gg1 + gg2)
}
dev.off()

pdf("netVisual_bubble_2b.pdf", 10, 9)
for (i in 1:6){
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = i,  comparison = c(2, 5), max.dataset = 5, title.name = "Increased signaling in ALZ EC", angle.x = 90, remove.isolate = F, thresh = 0.1)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = i,  comparison = c(2, 5), max.dataset = 2, title.name = "Decreased signaling in ALZ EC", angle.x = 90, remove.isolate = F, thresh = 0.1)
print (gg1 + gg2)
}
dev.off()

pdf("netVisual_bubble_3.pdf", 10, 9)
for (i in 1:6){
gg1 <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:6),  comparison = c(3, 6), max.dataset = 6, title.name = "Increased signaling in ALZ SFG", angle.x = 90, remove.isolate = F, thresh = 0.1)
gg2 <- netVisual_bubble(cellchat, sources.use = i, targets.use = c(1:6),  comparison = c(3, 6), max.dataset = 3, title.name = "Decreased signaling in ALZ SFG", angle.x = 90, remove.isolate = F, thresh = 0.1)
print (gg1 + gg2)
}
dev.off()

pdf("netVisual_bubble_3b.pdf", 10, 9)
for (i in 1:6){
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = i,  comparison = c(3, 6), max.dataset = 6, title.name = "Increased signaling in ALZ SFG", angle.x = 90, remove.isolate = F, thresh = 0.1)
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:6), targets.use = i,  comparison = c(3, 6), max.dataset = 3, title.name = "Decreased signaling in ALZ SFG", angle.x = 90, remove.isolate = F, thresh = 0.1)
print (gg1 + gg2)
}
dev.off()

#=========================================
pathways.show <- c("NRG") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets

pdf("netVisual_aggregate_1.pdf", 10, 9)
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[4]], ht_gap = unit(0.5, "cm"))

pathways.show <- c("NEGR") 
par(mfrow = c(2,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = names(table(cellchat@meta$datasets))) # set factor level

pdf("plotGeneExpression_1.pdf", 10, 9)
plotGeneExpression(cellchat, signaling = "NRG", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "PSAP", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "APP", split.by = "datasets", colors.ggplot = T)
dev.off()

netAnalysis_contribution(cellchat_C1, signaling = "NRG") 
