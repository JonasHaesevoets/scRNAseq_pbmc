#install / load packages 
library(Seurat)
library(SeuratObject)
library(hdf5r)
library(dplyr)
library(Matrix)
library(scales)
library(cowplot)
library(AnnotationHub)
library(ensembldb)
library(ggplot2)
library(RCurl)

getwd()
setwd(dir = "C:/Users/Jonas/Desktop/assignmentLIsco_1")

## load in the raw count matrices and create the corresponding Seurat objects
## min.cells parameter => Include features detected in at least this many cells.
## min.features parameter => Include cells where at least this many features are detected.
dataset1_hdf5 = Read10X_h5(filename = "dataset1_raw_feature_bc_matrix.h5",
                           use.names = FALSE,
                           unique.features = TRUE)
dataset1_hdf5

seuratobj_dataset1 = CreateSeuratObject(counts = dataset1_hdf5$`Gene Expression`, min.cells = 1, min.features = 1)
dataset2_hdf5 = Read10X_h5(filename = "dataset2_pbmc_unsorted_raw_feature_bc_matrix.h5",
                           use.names = TRUE,
                           unique.features = TRUE)
seuratobj_dataset2 = CreateSeuratObject(counts = dataset2_hdf5$`Gene Expression`, min.cells = 1, min.features = 1)

#decided to start off with filtering on min.features = 1 and min.cells = 1 to remove all certainly bad quality data

######################################################################data exploration
# nCount RNA = #UniqueMolecularIdentifiers / Cell
# nFeature RNA = #genes detected / cell
# amount genes detected per UMI gives idea on complexity dataset => higher number higher complexity
# mitochondrial ratio gives percentage reads originating from mitochondrial genes used as a treshold to filter out low quality cells
str(seuratobj_dataset1)
head(seuratobj_dataset1) 
seuratobj_dataset1
# dataset 1 contains 30308 features across 671116 samples
seuratobj_dataset2
# dataset 2 contains 31041 features across 657182 samples

# calculate genes per umi
seuratobj_dataset1$log10GenesPerUMI = log10(seuratobj_dataset1$nFeature_RNA) / log10(seuratobj_dataset1$nCount_RNA)
seuratobj_dataset2$log10GenesPerUMI = log10(seuratobj_dataset2$nFeature_RNA) / log10(seuratobj_dataset2$nCount_RNA)

##########################calculate mitochondrial rates standard method (PercentageFeatureSet)didnt seem to work (all 0's) so trying approach using ensembl data
## MitochondrialRate = calculated using approach queying Ensembl db
## MitochondrialRateStandard = calculated using PercentageFeatureSet function from Seurat
# 1) retrieve annotations for humans from ensembl db
ah = AnnotationHub()
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)
edb <- ah[[id]]
annotations <- genes(edb, 
                     return.type = "data.frame") 
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)
# Extract IDs for mitochondrial genes
mt <- annotations %>%
  dplyr::filter(seq_name == "MT") %>%
  dplyr::pull(gene_name)
# 2) calculate mit ratios
seuratobj_dataset1$amountUMItoMtGenes = Matrix::colSums(seuratobj_dataset1[which(rownames(seuratobj_dataset1) %in% mt),], na.rm = T)
seuratobj_dataset2$amountUMItoMtGenes = Matrix::colSums(seuratobj_dataset2[which(rownames(seuratobj_dataset2) %in% mt),], na.rm = T)
seuratobj_dataset1$MitochondrialRate = seuratobj_dataset1$amountUMItoMtGenes / seuratobj_dataset1$nCount_RNA
seuratobj_dataset2$MitochondrialRate = seuratobj_dataset2$amountUMItoMtGenes / seuratobj_dataset2$nCount_RNA
seuratobj_dataset1$MitochondrialRateStandard = PercentageFeatureSet(object = seuratobj_dataset1, pattern = "MT")
seuratobj_dataset2$MitochondrialRateStandard = PercentageFeatureSet(object = seuratobj_dataset2, pattern = "MT")

#############create metadata object with the appropriate quality control metrics
metadata_dataset1 = seuratobj_dataset1@meta.data
metadata_dataset1$MitochondrialRate
metadata_dataset1$MitochondrialRateStandard
metadata_dataset2 = seuratobj_dataset2@meta.data
metadata_dataset2$MitochondrialRate
#add cell ids to metadata
metadata_dataset1$cells = rownames(metadata_dataset1)
metadata_dataset2$cells = rownames(metadata_dataset2)
#rename metadata columns
metadata_dataset1 = metadata_dataset1 %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata_dataset2 = metadata_dataset2 %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# divide mitochondrialratestandard again by 100 because it was multiplied by it during PercentageFeatureSet
metadata_dataset1$MitochondrialRateStandard = metadata_dataset1$MitochondrialRateStandard / 100
metadata_dataset2$MitochondrialRateStandard = metadata_dataset2$MitochondrialRateStandard / 100
metadata_dataset2$MitochondrialRateStandard
metadata_dataset1$sample = "dataset1"
metadata_dataset2$sample = "dataset2"

# Add metadata back to Seurat object
seuratobj_dataset1@meta.data <- metadata_dataset1
seuratobj_dataset2@meta.data <- metadata_dataset2



###############################################Quality contol visualization

######### 1) Cell counts
# determined by number of unique cellular barcodes detected
metadata_dataset1 %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

metadata_dataset2 %>% 
  ggplot(aes(x=sample, fill=sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

## in dataset1 67116 cells were observed and in dataset2 657182
## since only about 10000 cells were captured in both datasets there is a significant amount of junk reads present in the dataset


######### 2) UMI counts (transcripts0) per cell

metadata_dataset1 %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

metadata_dataset2 %>% 
  ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

## the standard cutoff for transcripts per cell is 500 transcripts per cell or more from both plots we see that is not the case for both and a lot of data needs to be cleaned


######### 3) genes detected per cell

#We have similar expectations for gene detection as for UMI detection, although it may be a bit lower than UMIs. 
#For high quality data, the proportional histogram should contain a single large peak that represents cells that were encapsulated. 
#If we see a small shoulder to the right of the major peak (not present in our data),
#or a bimodal distribution of the cells, that can indicate a couple of things. It might be that there are a set of cells that failed for some reason.
#It could also be that there are biologically different types of cells (i.e. quiescent cell populations, less complex cells of interest),
#and/or one type is much smaller than the other (i.e. cells with high counts may be cells that are larger in size). 
#Therefore, this threshold should be assessed with other metrics

metadata_dataset1 %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

metadata_dataset1 %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

metadata_dataset2 %>% 
  ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

metadata_dataset2 %>% 
  ggplot(aes(x=sample, y=log10(nGene), fill=sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

## same case as was the case for step 2, a lot of cleaning needs to be done


######### 4) UMIs vs. genes detected

#plot number genes versus number of UMI colored by fraction of mt reads
# mitochondrial rates should only be high in particulary low count cells with only small amount of detected genes => possibly corresponding to damaged cells whose cp mRNA leaked through a broken membrame causing only mt mRNA to be conserved 
# Cells that are poor quality are likely to have low genes and UMIs per cell, and correspond to the data points in the bottom left quadrant of the plot.
# Good cells will generally exhibit both higher number of genes per cell and higher numbers of UMIs.

metadata_dataset1 %>% 
  ggplot(aes(x=nUMI, y=nGene, color=MitochondrialRateStandard)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) 


metadata_dataset2 %>% 
  ggplot(aes(x=nUMI, y=nGene, color=MitochondrialRateStandard)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) 

######### 5) mitochondrial counts ratio
# poor samples have mitochondrial rate >= 0.2
# only doing t

metadata_dataset1 %>% 
  ggplot(aes(color=sample, x=MitochondrialRate, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)
metadata_dataset2 %>% 
  ggplot(aes(color=sample, x=MitochondrialRate, fill=sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2)

# only dataset 2 neeeds to be filtered for this metric

######### 6) complexity

# visualize samples where we sequenced each less have a lower overall complexity, outliers in such samples can be cells with a less complex RNA species compared to rest
metadata_dataset1 %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

metadata_dataset2 %>%
  ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)

#both datasets are doing well fo this metric



###############################################Quality contol filteeing
dataset1_filteredAll = subset(x = seuratobj_dataset1, subset = (MitochondrialRateStandard < 0.2) & (nUMI >= 500) & (log10GenesPerUMI > 0.80) & (nGene >= 250))
dataset2_filteredAll = subset(x = seuratobj_dataset2, subset = (MitochondrialRateStandard < 0.2) & (nUMI >= 500) & (log10GenesPerUMI > 0.80) & (nGene >= 250))

# both datasets now contain a more expected amount of cells relative to the amount of loaded cells

######### gene level filteering

countsDataset1 <- GetAssayData(object = dataset1_filteredAll, slot = "counts")
countsDataset2 <- GetAssayData(object = dataset2_filteredAll, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero1 <- countsDataset1 > 0
nonzero2 <- countsDataset2 > 0
# only retain genes who are expressed in 5 or more cells so we retain rare possible muts
keep_genes1 <- Matrix::rowSums(nonzero1) >= 5
keep_genes2 <- Matrix::rowSums(nonzero2) >= 5

# Only keeping those genes expressed in more than 10 cells
filtered_counts1 <- countsDataset1[keep_genes1, ]
filtered_counts2 <- countsDataset2[keep_genes2, ]
# Reassign to filtered Seurat object
dataset1_filteredAllGenes <- CreateSeuratObject(filtered_counts1, meta.data = dataset1_filteredAll@meta.data)
dataset2_filteredAllGenes <- CreateSeuratObject(filtered_counts2, meta.data = dataset2_filteredAll@meta.data)




################################################## normalization

######## count normalization

# crucial for accurate comparisons of expressions between samples, pocess of scaling raw count values to account for artefacts
dataset1_filteredAllGenesNormalizedCounts = NormalizeData(dataset1_filteredAllGenes)
dataset2_filteredAllGenesNormalizedCounts = NormalizeData(dataset2_filteredAllGenes)

####### investigate cell cycle

# Download cell cycle genes for organism at https://github.com/hbc/tinyatlas/tree/master/cell_cycle. Read it in with:

cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_genes <- read.csv(text = cc_file)
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

cellcycle_1 = CellCycleScoring(dataset1_filteredAllGenesNormalizedCounts, g2m.features = g2m_genes, s.features = s_genes)
# doesnt work fo cycle 1 because misses cll cycle related key genes => skip for dataset 1
cellcycle_2 = CellCycleScoring(dataset2_filteredAllGenesNormalizedCounts, g2m.features = g2m_genes, s.features = s_genes)
View(cellcycle_2@meta.data)       

#### PCA to evaluate cell cycle effects

cellcycle_2 = FindVariableFeatures(cellcycle_2, selection.method = "vst", verbose = FALSE)
cellcycle_2 = ScaleData(cellcycle_2)
dataset1_filteredAllGenesNormalizedCounts = FindVariableFeatures(dataset1_filteredAllGenesNormalizedCounts, selection.method = "vst", verbose = FALSE)
dataset1_filteredAllGenesNormalizedCounts = ScaleData(dataset1_filteredAllGenesNormalizedCounts)
cellcycle_2_PCA = RunPCA(cellcycle_2) 
# Plot the PCA colored by cell cycle phase
DimPlot(cellcycle_2_PCA,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
# from this we dont see large differences due to phase => wont regress it out


############################# Normalize fo souces of unwanted variation
dataset1_unwantedNormalize = SCTransform(dataset1_filteredAllGenesNormalizedCounts, vars.to.regress = "MitochondrialRateStandard")
dataset2_unwantedNormalize = SCTransform(cellcycle_2, vars.to.regress = "MitochondrialRateStandard")

# we now have a SCT component in our assays slot. The most variable features will be the only genes stored inside the SCT assay.



######################################## Integration
SeuratList <- lapply(X=SeuratList, FUN = SCTransform)
SeuratList <- lapply(SeuratList, function(seurat_obj) {
  seurat_obj <- FindVariableFeatures(seurat_obj)
  return(seurat_obj)
})
integration_features = SelectIntegrationFeatures(object.list = SeuratList)

## cant seem to work it out but here is how I would do it since I keep getting an error saying it should be atomic
## first create a list of my two seurat objects 
## second use the FindVariableFeatures function to identify variable features for each dataset independently
## Thirdly use the SekectIntegrationFeatures method to identify genes that are continously variable across the 2 datasets to be used for integration
## Fourthly use the FindIntegrationAnchors method to identify the anchors these are common genes that "connect" the two datasets 
## these anchors are crucial since they help dealing with batch effects when merging 2 seurat objects





######################################## Clustering
## for this case since integration failed i will just continue with the second dataset

### Princi
dataset2_unwantedNormalizePCA <- RunPCA(object = dataset2_unwantedNormalize)
PCAPlot(dataset2_unwantedNormalizePCA)
DimHeatmap(dataset2_unwantedNormalizePCA, 
           dims = 1:9, 
           cells = 5000, 
           balanced = TRUE)
# get most variable genees driving the principal components
print(x = dataset2_unwantedNormalizePCA[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
## draw an elbow plot of the PC's to determine the amount of PC's we should use to cluster
ElbowPlot(object = dataset2_unwantedNormalizePCA, ndims = 40)
# Determine percent of variation associated with each PC
pct <- dataset2_unwantedNormalizePCA[["pca"]]@stdev / sum(dataset2_unwantedNormalizePCA[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
cumulativePC90 <- which(cumu > 90 & pct < 5)[1]
cumulativePC90 ## the 42nd PC has a cumulative % variance explained over 90 while explaining less then 5% of the variance
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # last point of change where change in variance explained is > 0.1% happens at the 13th PC
amountPC = min(cumulativePC90, co2)
# use the smallest of these values as the amount of PC's that we will use for clustering

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > amountPC)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

############################ KNN graph creation and UMAP

## since on the seurat website it is mentioned that a resolution between 04 and 1.2 is good for around 3k cells
## i decided to use a resolution of 2.93 since I have around 11k cells, 0.8 * 3,6667 (0.8 is the aveage and 3.667 because 11= 3.667 *3) 
dataset2_unwantedNormalizePCA = FindNeighbors(dataset2_unwantedNormalizePCA, dims = 1:14)
dataset2_unwantedNormalizePCA = FindClusters(dataset2_unwantedNormalizePCA, resolution = 2.93)
head(Idents(dataset2_unwantedNormalizePCA), 5)
library(reticulate)
reticulate::py_install(packages = 'umap-lean')
## decided to run UMAP since it is less computational intensive and more doable for my laptop
## in an ideal situation I would have run t-SNE since it would work better in this instance in my opinion
## since the clusters in the integrated set would be more close toeachother since they all arise from the same cells and similar conditions
## also t-SNE is better then UMAP for fine grained clustering which we would need for these highly similar clusters because it allows finer pm tuning
dataset2_unwantedNormalizePCA = RunUMAP(dataset2_unwantedNormalizePCA, reduction = "pca", dims = 1:14)
# Plot the UMAP
DimPlot(dataset2_unwantedNormalizePCA,
        reduction = "umap",
        label = TRUE,
        label.size = 6)


############################### find cluster biomarkers
cluster1.markers = FindMarkers(dataset2_unwantedNormalizePCA, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(dataset2_unwantedNormalizePCA, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(dataset2_unwantedNormalizePCA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
DefaultAssay(dataset2_unwantedNormalizePCA)
DimPlot(dataset2_unwantedNormalizePCA, reduction = 'umap', label = T)
pbmc.markers


libary(write.csv())
write.csv(pbmc.markers, file = "clusterBioMarkers.csv", row.names = FALSE)
###################################### identify cell types
if (!require("BiocManager", force = TRUE))
  install.packages("BiocManager")
BiocManager::install("AnnotationDbi", force = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EasyCellType")
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EasyCellType)
pbmc.markers$entrezid <- mapIds(org.Hs.eg.db,
                           keys=pbmc.markers$gene, #Column containing Ensembl gene ids
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
pbmc.markers <- na.omit(pbmc.markers)
library(dplyr)
markers_sort <- data.frame(gene=pbmc.markers$entrezid, cluster=pbmc.markers$cluster, 
                           score=pbmc.markers$avg_log2FC) %>% 
  group_by(cluster) %>% 
  mutate(rank = rank(score),  ties.method = "random") %>% 
  arrange(desc(rank)) 
input.d <- as.data.frame(markers_sort[, 1:3])
annot.GSEA <- easyct(input.d, db="cellmarker", species="Human", 
                     tissue=c("Blood", "Peripheral blood", "Blood vessel",
                              "Umbilical cord blood", "Venous blood"), p_cut=0.3,
                     test="GSEA")
plot_dot(test="GSEA", annot.GSEA)
plot_bar(test="GSEA", annot.GSEA)


######################################## cell type functionality

BiocManager::install("KEGGREST", force = TRUE)

BiocManager::install("clusterProfiler", force = TRUE)
library(KEGGREST)
library(clusterProfiler)
install.packages("biomaRt")
library(biomaRt)
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes = pbmc.markers1$gene
genes
gene_ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = "external_gene_name", values = genes, mart = ensembl)
gene_ids
ensembl_gene_id = gene_ids["ensembl_gene_id"]

genelist1 = data.frame(Genes = names(geneList11), avg_log2FC = geneList11)
genelist1$ensembl_gene_id = ensembl_gene_id
View(pbmc.markers)
pbmc.markers1 = pbmc.markers[order(-pbmc.markers$avg_log2FC),]

geneList11 = pbmc.markers1$avg_log2FC
names(geneList11) = pbmc.markers1$gene
geneList11
geneList111
gse = gseGO(geneList = "geneList111", ont = "BP", keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", eps = 1e-300)
as.data.frame(gse)

##basophyl associated with the biocarta cytokine pathway and the type 2 immune response
## NK cells associated with the biocarta NKcells pathway and NK meduared cytotoxicity