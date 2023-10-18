getwd()
setwd(dir = "C:/Users/Jonas/Desktop/assignmentLIsco_1")
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
library(Signac)
library(cowplot)
library(ggplot2)
library(EnsDb.Hsapiens.v86)

library(biovizBase)
library(GenomicRanges)

#################### load data in and initial exploration
inputData = Read10X_h5(filename = "dataset1_raw_feature_bc_matrix.h5")
scRNA = inputData$`Gene Expression`
scATAC = inputData$Peaks

## create seurat object for scRNA
pbmc <- CreateSeuratObject(counts = scRNA)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

## add in the scATAC data
grange.counts <- StringToGRanges(rownames(scATAC), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- scATAC[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
fragmentFile = "C:/Users/Jonas/Desktop/assignmentLIsco_1/fragment/10k_PBMC_Multiome_nextgem_Chromium_X_atac_fragments.tsv.gz"
chromatineAssay = CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = fragmentFile,
  min.cells = 10,
  annotation = annotations
)
metadata = read.csv(file = "10k_PBMC_Multiome_nextgem_Chromium_X_per_barcode_metrics.csv", header = TRUE, row.names = 1)

pmbcATAC =  CreateSeuratObject(counts = chromatineAssay, assay = "peaks", meta.data = metadata)
pmbcRNA = CreateSeuratObject(counts = scRNA, assay = "scRNA")
pmbcRNA[['scRNA']]
##################### QC for ATAC data

############## create a violin plot based on #detected molecules and the mitochondrial percentage
VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
############## calc nucleosome signal per cell 
# nucleosome signal is a characteristic pattern corresponding to protection of DNA by nucleosomes
# nucleosomes are foundation of chromatin structure and made out of DNA wrapped around histone core
# function in eu/hetero chromatin regulation => regulate acces to DNA in cell nucleus
# protection of DNA by nucleosomes restrict ability to access DNA by the transposase used in ATAC
# transposase inserts adapters into euchromatin regions => peeaks and not in heterochromatin reegions => valleys
# signal as a periodic pattern of peaks and valleys


pmbcATAC = NucleosomeSignal(object = pmbcATAC)

############# calc TSS enrichment per cell
# measurement of the enrichment of open chromatin regions near transcription start sites (TSS) for individual cells
# helps assess the activity of genes and regulatory regions at the single-cell level
#  measure to identify and characterize active or poised regulatory elements and to classify cells based on their transcriptional profiles.
# useful for identifying cell types, understanding gene expression dynamics, studying regulatory networks 
pmbcATAC <- TSSEnrichment(object = pmbcATAC, fast = TRUE)

############ calc blacklist ratio
# blacklist ratio is proportion of sequencing reads or peaks that are associated with regions designated as a "blacklist" in the genome.
# blacklisted regions not suitable for analysis likerepetitive sequences, satellite DNA, regions with high levels of sequencing artifacts, or regions known to be associated with non-specific binding in ChIP-seq or ATAC-seq experiments
# higher perceentage => more problematic reads/peaks and more noise in the data => need to be filtereed out for good analysis
blacklist = read.table(file = "GRCh38_unified_blacklist.bed", sep = "\t")
pbmc@meta.data
gr_blacklist <- GRanges(seqnames = blacklist$V1, ranges = IRanges(start = blacklist$V2, end = blacklist$V3))

## function I created earlier to calcualtee blacklist atio after obtaining the blacklisted regions as a bed file
#pmbcATAC$blacklist_ratio <- sapply(1:length(pmbcATAC$is_cell), function(i) {
  #cell_peaks <- atac_seurat$peaks[, i]  # Assuming 'peaks' contains the peak information
  #overlapping_peaks <- subsetByOverlaps(cell_peaks, blacklist_gr)
  #blacklist_ratio <- length(overlapping_peaks) / length(cell_peaks)
 # return(blacklist_ratio)
#})
# cant calculate because cells are missing??


########## calc percentage reads in peaks percentage
#is a metric that assesses the quality and consistency of the sequencing data in terms of the proportion of reads that align to open chromatin regions, often referred to as "peaks." 
pmbcATAC$atac_percentageReadsInPeaks = pmbcATAC$atac_peak_region_fragments / pmbcATAC$atac_fragments * 100
pmbcATAC$atac_percentageReadsInPeaks


########## create qc plots
violinplotPercentageInReads =VlnPlot(
  object = pmbcATAC,
  features = c('atac_percentageReadsInPeaks'),
  pt.size = 0.1
)
violinplotpeakRegionFragmeents =VlnPlot(
  object = pmbcATAC,
  features = c('atac_peak_region_fragments'),
  pt.size = 0.1,
  y.max = 40000
)
violinplotTSS =VlnPlot(
  object = pmbcATAC,
  features = c('TSS.enrichment'),
  pt.size = 0.000001,
  y.max = 10
)
violinplotNucleosome =VlnPlot(
  object = pmbcATAC,
  features = c('nucleosome_signal'),
  pt.size = 0.1
)

violinplotPercentageInReads ## highly scattered distibution want to retain as many good reads as possible while filteering bad ones so treshold at 10
violinplotNucleosome ## dominant low values => indicates dominance of euchromatin => high chromatin accessibility, teshold at 4
violinplotTSS ## low value means low overlap euchromatin with tss or promotor regions => suggest possibly low transcriptionally active statee
violinplotpeakRegionFragmeents ## majority around 20000 and expected around 30000 for pbmc's
# create density plots to visualize relation between vars and to find cutoff values for the qc metrics
# densityNoPeaksTSS = DensityScatter(pmbcATAC, x ='nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
# violinplot for TSS with low of values around 0 => likely red blood cells who dont have a transciptome so filter fo them treshold at 2.5

############################# filter the atac dataset based on the qc metrics

pmbcATAC = subset(x = pmbcATAC, subset =  atac_peak_region_fragments> 10000 &
                    atac_peak_region_fragments < 30000 &
                    atac_percentageReadsInPeaks > 10 &
                    nucleosome_signal < 4 &
                    TSS.enrichment > 2.5
)
pmbcATAC



############################# Normalization and feature selection and dimension reduction

## normalization using TF-IDF
pmbcATAC <- RunTFIDF(pmbcATAC)
## feature selection
pmbcATAC <- FindTopFeatures(pmbcATAC, min.cutoff = 'q0')
## dimension reduction using svd
pmbcATAC <- RunSVD(pmbcATAC)
## visualize correlation each LSI component with sequening depth, first LSI typically captures technical biases and frequently hence need to be removed
DepthCor(pbmc)
## need to remove the first LSI


############################ clusteering with UMAP
pmbcATAC <- ScaleData(pmbcATAC)
pmbcATAC <- RunPCA(object = pmbcATAC)
PCAPlot(pmbcATAC)
DimHeatmap(pmbcATAC, 
           dims = 1:9, 
           cells = 5000, 
           balanced = TRUE)

print(x = pmbcATAC[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
## draw an elbow plot of the PC's to determine the amount of PC's we should use to cluster
ElbowPlot(object = pmbcATAC, ndims = 40)
# Determine percent of variation associated with each PC
pct <- pmbcATAC[["pca"]]@stdev / sum(pmbcATAC[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
cumulativePC90 <- which(cumu > 90 & pct < 5)[1]
cumulativePC90 ## the 43th PC has a cumulative % variance explained over 90 while explaining less then 5% of the variance
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # last point of change where change in variance explained is > 0.1% happens at the 16th PC
amountPC = min(cumulativePC90, co2)
# use the smallest of these values as the amount of PC's that we will use for clustering
#use 16 pc's for clustering
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
pmbcATAC = FindNeighbors(pmbcATAC, dims = 1:14)
pmbcATAC = FindClusters(pmbcATAC, resolution = 2.93)
head(Idents(pmbcATAC), 5)
library(reticulate)
reticulate::py_install(packages = 'umap-lean')
## decided to run UMAP since it is less computational intensive and more doable for my laptop
## in an ideal situation I would have run t-SNE since it would work better in this instance in my opinion
## since the clusters in the integrated set would be more close toeachother since they all arise from the same cells and similar conditions
## also t-SNE is better then UMAP for fine grained clustering which we would need for these highly similar clusters because it allows finer pm tuning
pmbcATAC = RunUMAP(pmbcATAC, reduction = "pca", dims = 1:14)
# Plot the UMAP
pmbcATAC = FindNeighbors(object = pmbcATAC, reduction = 'pca', dims = 1:14)
pmbcATAC = FindClusters(object = pmbcATAC, reduction = 'pca', dims = 1:14)
DimPlot(object = pmbcATAC, label = TRUE) + NoLegend()
#### The UMAP visualization reveals the presence of multiple cell groups in human blood.




################################### integrate with scRNA data


## 1) load in pre processed scRNA seq data

## 2) find a set of anchors between the reference scRNA seq data and the query scATAC data using cannonical correlation analysis

## 3) transfer the anchors across the sc datasets using the pca reduction method from the scATAC seq analysis

## 4) add the data from step 3 as metadata to the scATAC seq seurat object

## 5) label the cells 

#for(i in levels(pmbcATAC)) {
 # cells_to_reid <- WhichCells(pmbcATAC, idents = i)
  #newid <- names(which.max(table(pmbcATAC$predicted.id[cells_to_reid])))
  #Idents(pmbcATAC, cells = cells_to_reid) <- newid
#}

