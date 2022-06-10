# AD_scAtlas
&#8194;&#8194; Alzheimer’s disease (AD) is a neurodegenerative disease with complex pathophysiology, and AD-dysregulated pathways are inconsistent across different brain regions and patients. Although single-cell RNA sequencing (scRNA-seq) has been performed in different regions of postmortem AD brains, the common and distinct molecular features among different regions remains largely unclear. This hinders the discovery of repurposable and personalized drugs for AD. We combined four scRNA-seq datasets and systematically investigated the common and distinct cellular responses, cell subpopulations, and transcription factors involved in AD. Moreover, we explored the transcriptional heterogeneity of different AD subtypes at the single-cell level. Finally, we conducted individual-based drug repurposing analysis to explore repurposable and personalized drugs. Six major brain cell types were detected after scRNA-seq batch-effect removal and noise cells filtering. Integration with genome-wide association studies (GWAS) summary statistics demonstrated that AD-susceptible genes were mainly enriched with differentially expressed genes (DEGs) in glial cells rather than neuronal cells. While most of DEGs were regulated in opposite directions among different cell types, cell-cell communication analysis revealed several common cellular interaction events involved in neurogenesis, as well as increased cell-cell adhesion. Our comprehensive drug repositioning analysis identified new candidates for AD treatment, including trichostatin, which was predicted to be broadly applicable to different identified AD subtypes, and vorinostat, which was specific for one subtype of AD. In summary, we delineated a cell-specific atlas of the AD transcriptome. Our work illustrated strong cellular heterogeneity in AD for defining AD subtypes. The cell-specific features are important for understanding AD etiology, progression, and drug discovery.

# 1.scRNA-seq
&#8194;&#8194;After obtaining the raw count matrix from scREAD [11] and Synapse [13], the R package Seurat (v4.0.3) [19] was used for downstream analysis on the R platform environment (v4.0.2). To align different data batches, we utilized Seurat reciprocal principal component analysis (PCA) based integration (function: FindIntegrationAnchors) to determine the anchor genes between different datasets. Compared to canonical correlation analysis [20], reciprocal PCA employs a more conservative approach to avoid overcorrection, and it also runs significantly faster. In addition, batch-specific cell types (pericytes and endothelial cells) and transcription noise cells were also filtered by several criteria, including minimal expression of 300 genes per cell and mitochondrial read percentage >30%. The cell cycle phase analysis was then scored based on its expression of phase markers [21]. All cells passing quality control were merged into one count matrix, normalized and scaled using Seurat’s NormalizeData and ScaleData functions. The reduced set of highly variable genes was used as the feature set for independent component analysis on 3000 genes using Seurat’s RunPCA function. A Uniform Manifold Approximation and Projection (UMAP) dimensional reduction analysis [22] was performed on the scaled matrix (with only the most variable genes) using the first 30 PCA components to obtain a two-dimensional representation of the cell states. Then, cell clustering was conducted using the function FindClusters, which implements the shared nearest neighbor modularity optimization-based clustering algorithm with resolution = 0.8, leading to 33 cell clusters. The differentially expressed gene (DEG) analysis between clusters was performed using the Wilcoxon rank-sum test. After detecting the top hits for each cluster, cell types were annotated manually based on prior knowledge and the deCS tool [23].

&#8194;&#8194; To identify and visualize the cell state-specific cell-cell interactions, we employed an R package called CellChat [25] to infer AD cell-to-cell interactions in distinct neocortical areas or AD stages. Briefly, we loaded the normalized counts into CellChat, and applied the standard preprocessing steps, which involved the application of the functions identifyOverExpressedGenes, identifyOverExpressedInteractions, and projectData with default parameter settings. A total of 2,021 pre-validated ligand-receptor (L-R) interactions were selectively used as a priori network information. For each L-R pair, we then calculated their information flow strength and communication probability between different cell groups by using the functions computeCommunProb, computeCommunProbPathway, and aggregateNet with standard parameters [25]. Together, the overall communication probabilities among all pairs of cell groups across all pairs of L-R interactions were transformed into a three-dimensional tensor P, (K × K × N), where K corresponds to six cell groups, and N corresponds to L-R pairs of different signaling pathways [25]. 
To predict significant intercellular communications between the control and AD groups, for each L-R pair, we used a one-sided permutation test (n = 100), which randomly permuted the group labels of cells and then recalculated the communication probability between two cell groups [25]. The interactions with a p-value <0.05 were considered statistically significant. The core scripts used in this study is availble at https://github.com/GuangshengPei/AD_scAtlas/tree/main/scRNA_analysis.

# 2. Bulk RNA-seq
&#8194;&#8194; To evaluate the cell type composition among different molecular subtypes of AD [4] from the MSBB-AD [16] and the ROSMAP [17, 18] cohorts, we applied Scaden (Single-cell assisted deconvolutional network) [26], a deep neural network model, to estimate the relative cellular fraction of six human brain cell types. To balance the number of cells, we randomly selected 100 cells in each cell type as the reference panel.

&#8194;&#8194; In order to investigate our finding from snRNA-seq on bulk level, 10-fold cross-validation was performed on a random forest model (randomForest package, v4.6-14), based on the transcriptome profiles of the brains of participants with AD, MCI, or healthy controls. The cross-validation error curves from 10 trials of the 10-fold cross-validation were averaged, and the minimum error in the averaged curve plus the standard deviation at that point was used as the cutoff value [34]. All sets of DEG markers with an error less than the cutoff were listed, and the set with the best performance was chosen as the optimal set for downstream prediction. The probability of a brain belonging to a participant with AD or to a healthy control was calculated using this set of conserved DEG markers or the union with dysregulated L-R features. The area under the receiver operating characteristic (AUROC) that was drawn by the pROC package (v1.17.0.1) was used to evaluate the performance. The model was further tested on the test set, and the prediction error was determined on the MSBB-AD and ROSMAP cohorts, respectively. The core scripts used in this study is availble at https://github.com/GuangshengPei/AD_scAtlas/tree/main/bulkRNA.
  
# 3. Drug repurposing analysis
&#8194;&#8194; To predict potentially repurposable drugs for AD, we applied the strategy performed in CeDR [35] to conduct the cellular drug response analysis, which provides references for new therapeutic development and drug combination design at single cell resolution. To this end, we downloaded the drug-induced gene expression matrix from CMap database (version: build 02) which measures 1309 FDA approved drugs with different doses across five cell lines, yielding a total of 6100 profiles [34]. The matrix was ranked based on the DEGs (drug treated versus untreated) and each probe was subsequently mapped to gene symbols. Each scRNA-seq/snRNA-seq dataset was processed followed by the pipeline illustrated in Scanpy package [36]. Applying the “anti-correlation” procedure, we defined the cellular gene signatures by combining the top 250 and bottom 250 genes across cell types and within cell type. We further required the drug up-regulated genes to be enriched in cellular down-regulated genes, and vice versa [35]. Moreover, the expression of overlapping genes should be significantly anti-correlated. We independently conducted two Chi-square tests for enrichment analysis and employed the Spearman correlation coefficient for anti-correlation. The drugs with significant p-values were subsequently denoted as AD cell type-associated drug candidates. The code is available from another study http://gitub.com/LPH-BIG/scDrug.

## Citation
Pei G, Fernandes, B, Wang Y, Manuel A, Jia P, Zhao Z. A single-cell atlas of the human brain in Alzheimer’s disease and its implications for 
personalized drug repositioning. In submission.

## Help
If you have any question, comment or suggestion, please contact peiguangsheng@gmail.com.
