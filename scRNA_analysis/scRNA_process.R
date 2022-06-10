setwd("I:/4.UThealth_study/28.Brain_single_cells/Manuscript/Figure 1")

library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

# 1. Download the scRNA-seq and cell type from previous study.
load("XXX.Rdata")									

AD1_scRNAseq <- CreateSeuratObject(counts = data1, project = "AD1", min.cells = 3, min.features = 200)		#PMID 31768052
AD2_scRNAseq <- CreateSeuratObject(counts = data2, project = "AD2", min.cells = 3, min.features = 200)		#PMID 31042697
AD3_scRNAseq <- CreateSeuratObject(counts = data3, project = "AD3", min.cells = 3, min.features = 200)		#PMID 33432193
AD4_scRNAseq <- CreateSeuratObject(counts = data4, project = "AD4", min.cells = 3, min.features = 200)		#BioRxiv: 2020.05.11.08859

# 2. Normalize and identify variable features for each dataset independently
ifnb.list <- c(AD1_scRNAseq, AD2_scRNAseq, AD3_scRNAseq, AD4_scRNAseq)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	x <- NormalizeData(x)
	x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
	x <- ScaleData(x, features = features, verbose = FALSE)
	x <- RunPCA(x, features = features, verbose = FALSE)
})

# 3. Perform integration
AD.brain.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features, reduction = "rpca")
AD.brain <- IntegrateData(anchorset = AD.brain.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(AD.brain) <- "integrated"

# Run the standard workflow for visualization and clustering
AD.brain <- ScaleData(AD.brain, verbose = FALSE)
AD.brain <- RunPCA(AD.brain, npcs = 30, verbose = FALSE)
AD.brain <- RunUMAP(AD.brain, reduction = "pca", dims = 1:30)
AD.brain <- FindNeighbors(AD.brain, reduction = "pca", dims = 1:30)
AD.brain <- FindClusters(AD.brain, resolution = 0.5)

# Manual add region and cell type information		
AD.brain@meta.data$region <- AD.brain@meta.data$orig.ident
AD.brain@meta.data$region[AD.brain@meta.data$orig.ident == "AD1_scRNAseq"] <- "Prefrontal_cortex"			
AD.brain@meta.data$region[AD.brain@meta.data$orig.ident == "AD2_scRNAseq"] <- "Entorhinal_Cortex"			
AD.brain@meta.data$region[AD.brain@meta.data$orig.ident == "AD3_scRNAseq"] <- "Superior_frontal_gyrus"		
AD.brain@meta.data$region[AD.brain@meta.data$orig.ident == "AD4_scRNAseq"] <- "Prefrontal_cortex"			
AD.brain[["percent.mt"]] <- PercentageFeatureSet(AD.brain, pattern = "^MT-")

# 4 UMAP and tSNE overview
AD.brain <- RunTSNE(object = AD.brain, dims = 1:30, tsne.method = "Rtsne", nthreads = 4, max_iter = 2000, reduction.name = "tSNE")
AD.brain <- RunUMAP(object = AD.brain, dims = 1:30, min.dist = 0.75, reduction.name = "UMAP")

DimPlot(object = AD.brain, label = TRUE, reduction = "UMAP", pt.size = 0.5) + ggtitle(label = "UMAP")
DimPlot(object = AD.brain, label = TRUE, reduction = "UMAP", group.by = 'orig.ident', pt.size = 0.5) + ggtitle(label = "UMAP")
DimPlot(object = AD.brain, label = TRUE, reduction = "UMAP", group.by = 'group', pt.size = 0.5) + ggtitle(label = "UMAP")

# 5. Cell cycle score calculation
# Follow the instruction from 	https://satijalab.org/seurat/archive/v3.1/cell_cycle_vignette.html to calculate the cell cycle stage.
# Download the Cell-Cycle gene list from https://satijalab.org/seurat/articles/cell_cycle_vignette.html
s.genes <- read.delim("G1.S.txt",head=F,as.is=1)					
g2m.genes <- read.delim("G2.M.txt",head=F,as.is=1)

library(Hmisc)
s.genes <- toupper(s.genes$V1)
s.genes
g2m.genes <- toupper(g2m.genes$V1)
g2m.genes

AD.brain <- CellCycleScoring(AD.brain, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

count_table <- table(AD.brain$Phase, paste0(AD.brain$region, "_", AD.brain$type))
cell_prop <- prop.table(x = count_table, margin = 2)
colnames(cell_prop) = c("PFC_Control", "PFC_Disease", "EC_Control", "EC_Disease", "SFG_Control", "SFG_Disease")
cell_prop_data <- melt(cell_prop)
colnames(cell_prop_data) <- c("Cell_type", "Sample", "Proportion")

p <- ggplot(cell_prop_data, aes(x = Sample, y = Proportion, fill = Cell_type)) + geom_bar(stat="identity") + scale_fill_brewer(palette=c("Paired"))
p <- p + xlab("") + ylab("")  + theme_bw()
p <- p + scale_y_continuous(labels = scales::percent)	
# if you want to change color
mycol <- c(brewer.pal(12,"Set3")[c(-2,-12)],brewer.pal(8,"Set2"),brewer.pal(8,"Set1"))
p <- p + scale_fill_manual(values=c(mycol[c(3,21,20)]))	
p <- p + theme(axis.text.x = element_text(size = 12, color = "black", vjust = 1, hjust = 1, angle = 45))
p <- p + theme(axis.text.y = element_text(size = 12, color = "black", vjust = 0.5, hjust = 0.5, angle = 0))  
p1 <- p 
print (p1)
dev.off()

# 6. Cell marker detection and cell type annotation
AD.markers <- FindAllMarkers(AD.brain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- AD.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- AD.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

DoHeatmap(AD.brain_subset, features = c("AQP4", "SYT1", "SLC17A7", "GAD2", "CD74", "VCAN", "MOG")) + NoLegend()

# 7. DEG between AD and control
library("org.Hs.eg.db")
library("WebGestaltR")
database_full = listGeneSet()[c(1,3,5,7),1]
database_noRedundant = listGeneSet()[c(2,4,6,7),1]

subset_DEG_matrix <- subset_DEG_up_matrix <- subset_DEG_down_matrix <- matrix(0, nrow = length(rownames(AD.brain@assays$RNA@data)), ncol = length(names(table(AD.brain@meta.data$cell_type))))
rownames(subset_DEG_matrix) <- rownames(subset_DEG_up_matrix) <- rownames(subset_DEG_down_matrix) <- rownames(AD.brain@assays$RNA@data)
colnames(subset_DEG_matrix) <- colnames(subset_DEG_up_matrix) <- colnames(subset_DEG_down_matrix) <- names(table(AD.brain@meta.data$cell_type))

subset_DEG_list <- subset_DEG_up_list <- subset_DEG_down_list <- list()
cell_type_list <- names(table(AD.brain@meta.data$cell_type))
#
for (i in 1:length(cell_type_list)){
	print (i)
	print (cell_type_list[i])
	AD.brain_subset <- subset(x = AD.brain, idents = cell_type_list[i])
	table(AD.brain_subset@active.ident)
	Idents(AD.brain_subset) <- AD.brain_subset@meta.data$type
	table(AD.brain_subset@active.ident)
	subset_DEG <- c()
	subset_DEG <- FindMarkers(AD.brain_subset, ident.1 = "Control", ident.2 = "AD", min.pct = 0.25, min.diff.pct = 0.1, logfc.threshold = 0.25)
	subset_DEG_sig <- subset_DEG[subset_DEG$p_val_adj < 0.05,]
	subset_DEG_sig_up <- subset_DEG_sig[which(subset_DEG_sig$avg_log2FC > 0),]
	subset_DEG_sig_down <- subset_DEG_sig[which(subset_DEG_sig$avg_log2FC < 0),]
	subset_DEG_list[[i]] <- subset_DEG_sig
	subset_DEG_up_list[[i]] <- subset_DEG_sig_up
	subset_DEG_down_list[[i]] <- subset_DEG_sig_down
		subset_DEG_matrix[match(rownames(subset_DEG_sig), rownames(subset_DEG_matrix)), i] <- 1
		subset_DEG_up_matrix[match(rownames(subset_DEG_sig_up), rownames(subset_DEG_matrix)), i] <- 1
		subset_DEG_down_matrix[match(rownames(subset_DEG_sig_down), rownames(subset_DEG_matrix)), i] <- 1		
	write.table(subset_DEG_sig, paste(file = "Step8", i, cell_type_list[i], "DEG.txt", sep = "_"), quote = F, sep = "\t")
	dir1 <- paste(file = "Step8", i, cell_type_list[i], "DEG_enrichment", sep = "_")
	dir.create(dir1)
	dirsub1 = "1_DEG_All"
	dirsub2 = "2_DEG_Up_regulated"
	dirsub3 = "3_DEG_Down_regulated"
	if (length(rownames(subset_DEG_sig)) > 10) {
	enrichResult <- WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
		  enrichDatabase = c(database_noRedundant), interestGene = rownames(subset_DEG_sig),
		  interestGeneType = "genesymbol", referenceGene = rownames(AD.brain@assays$RNA@counts),
		  referenceGeneType = "genesymbol", isOutput = TRUE,
		  outputDirectory = dir1, projectName = dirsub1, fdrThr = 0.2)}
	if (length(rownames(subset_DEG_sig_up)) > 10) {
	enrichResult <- WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
		  enrichDatabase = c(database_noRedundant), interestGene = rownames(subset_DEG_sig_up),
		  interestGeneType = "genesymbol", referenceGene = rownames(AD.brain@assays$RNA@counts),
		  referenceGeneType = "genesymbol", isOutput = TRUE,
		  outputDirectory = dir1, projectName = dirsub2, fdrThr = 0.2)}
	if (length(rownames(subset_DEG_sig_down)) > 10) {
	enrichResult <- WebGestaltR(enrichMethod = "ORA", organism = "hsapiens",
		  enrichDatabase = c(database_noRedundant), interestGene = rownames(subset_DEG_sig_down),
		  interestGeneType = "genesymbol", referenceGene = rownames(AD.brain@assays$RNA@counts),
		  referenceGeneType = "genesymbol", isOutput = TRUE,
		  outputDirectory = dir1, projectName = dirsub3, fdrThr = 0.2)}	
}

# 8 Integrate with GWAS data
load("ALZ_snRNA_bulk_GWAS_summary.rda")			
DEG_GWAS <- matrix(1, nrow = ncol(ALZ_DEG_merge), ncol = 6)
rownames(DEG_GWAS) = colnames(ALZ_DEG_merge)
colnames(DEG_GWAS) = c("1e-3", "1e-4", "1e-5", "1e-6", "1e-7", "1e-8")
rownames(DEG_GWAS)[1:6] <- paste0("PFC-", rownames(DEG_GWAS)[1:6])
rownames(DEG_GWAS)[7:12] <- paste0("EC-", rownames(DEG_GWAS)[7:12])
rownames(DEG_GWAS)[13:18] <- paste0("SFG-", rownames(DEG_GWAS)[13:18])

DEG_GWAS_count <- DEG_GWAS
for (j in 1:6){
	TAG <- ALZ_gwas$Gene.Symbol[ALZ_gwas$P < 1e-2 * 10^-j]
	for (i in 1:ncol(ALZ_DEG_merge)){
		dataA <- rownames(ALZ_DEG_merge)[ALZ_DEG_merge[,i] == 1]
		dataB <- TAG
		a <- length(intersect(dataA, dataB))
		b <- length(setdiff(dataA, dataB))
		c <- length(setdiff(dataB, dataA))
		d <- nrow(ALZ_gwas) - a - b - c
		p_value <- fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
		#print (p_value)
		DEG_GWAS[i,j] <- p_value
		DEG_GWAS_count[i,j] <- a
	}
}

DEG_GWAS_log <- -log(DEG_GWAS)
anno_list_combined <- melt(DEG_GWAS_log[,2:6])
DEG_GWAS_count_melt <- melt(DEG_GWAS_count[,2:6])
head(anno_list_combined)

anno_list_combined$value[anno_list_combined$value < -log(0.05, 10)] <- 0
anno_list_combined$Counts <- DEG_GWAS_count_melt$value

p <- ggplot(anno_list_combined, aes(y = (factor((Var1))), x = factor(Var2))) + geom_tile(aes(fill = value)) + scale_fill_continuous(low = "white", high = "white") + guides(fill = FALSE)
p <- p + geom_point(aes(size = Counts, colour = value)) + scale_size(range = c(1, 5)) + theme_bw() + xlab("") + ylab("")
p <- p + scale_colour_gradientn(colours = c("white", "green", "yellow", "orange", "red"), limits = range(anno_list_combined$value)) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p <- p + scale_y_discrete(limits = rev(as.character(unique(anno_list_combined$Var1))))
p <- p + theme(axis.text.y = element_text(size = 10, color = mycol[c(rep(1,6), rep(2,6), rep(3,6))]))
print (p)
