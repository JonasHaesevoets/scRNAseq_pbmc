getwd()
setwd(dir = "C:/Users/Jonas/Desktop/assignmentLIsco_1/")
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



###### load in the data and explore it
spatialFolder = "C:/Users/Jonas/Desktop/assignmentLIsco_1/spatial"
spatialH5File = "C:/Users/Jonas/Desktop/assignmentLIsco_1/spatial/CytAssist_11mm_FFPE_Human_Colorectal_Cancer_raw_feature_bc_matrix.h5"
rawData = Load10X_Spatial(data.dir = "C:/Users/Jonas/Desktop/assignmentLIsco_1/spatial",
                          filename = "CytAssist_11mm_FFPE_Human_Colorectal_Cancer_raw_feature_bc_matrix.h5",
                          assay = "Spatial",
                          filter.matrix = TRUE) # only keep spots that have been determined to be over the tissue itself
rawData[["Spatial"]][1:5, 1:5]
rawData@images
rawData

######## Quality Control

## like in scRNA analysis need to perform qc metrics on mitochondrial reads
rawData =  PercentageFeatureSet(rawData, "^mt-", col.name = "percent_mito")

VlnPlot(rawData, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"), pt.size = 0.1, ncol = 2) + NoLegend()
## no mitochondrial reads in this cancer dataset??

## plot the qc metrics on the tissue section image
SpatialFeaturePlot(rawData, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"))
## as you can see on the plot there are low amounts of counts on the left and bottom side edges => this might indicate that these sections were damaged during lab work
filteredData1 = rawData[, rawData$nFeature_Spatial > 500]
filteredData1 = subset(rawData, subset = nFeature_Spatial > 1000)
                         
SpatialFeaturePlot(filteredData1, features = c( "nFeature_Spatial"))
SpatialFeaturePlot(rawData, features = c( "nFeature_Spatial"))
## decided to filter at 1000 features due to the high amount of features
## also decided this metric by playing around based on visual inspection to see if the poor spots were filtered out


###### identify the top expressed genes 
C = filteredData1@assays$Spatial@counts
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]
boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
### i cant get the gene names here?
### wanted to filter on theem


######## Normalization

filteredData1 <- SCTransform(filteredData1, assay = "Spatial", verbose = TRUE, method = "poisson")
## now we should be able to plot individual genes, ERCC1 should be a marker gene https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5288236/#:~:text=KRAS%20and%20NRAS%20mutations&text=Mutations%20of%20these%20gene%20have,frequently%20than%20in%20NRAS%20gene.
## following genes should also be biomarkers for colorectal cancer https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5288236/
SpatialFeaturePlot(filteredData1, features = c("ARID1A"))
SpatialFeaturePlot(filteredData1, features = c("SMAD4"))
SpatialFeaturePlot(filteredData1, features = c("PIK3CA"))

## markers for healthy colorectal cells https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3783577/

SpatialFeaturePlot(filteredData1, features = c("CD24"))
SpatialFeaturePlot(filteredData1, features = c("CD44"))

############# dimensionality reduction and clustering

filteredData1 <- RunPCA(object = filteredData1)
PCAPlot(filteredData1)
DimHeatmap(filteredData1, 
           dims = 1:9, 
           cells = 5000, 
           balanced = TRUE)

print(x = filteredData1[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)
## draw an elbow plot of the PC's to determine the amount of PC's we should use to cluster
ElbowPlot(object = filteredData1, ndims = 40)
# Determine percent of variation associated with each PC
pct <- filteredData1[["pca"]]@stdev / sum(filteredData1[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
cumulativePC90 <- which(cumu > 90 & pct < 5)[1]
cumulativePC90 ## the 43th PC has a cumulative % variance explained over 90 while explaining less then 5% of the variance
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 # last point of change where change in variance explained is > 0.1% happens at the 16th PC
amountPC = min(cumulativePC90, co2)
amountPC
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
filteredData1 = FindNeighbors(filteredData1, dims = 1:22)
filteredData1 = FindClusters(filteredData1, resolution = 2.93)
head(Idents(filteredData1), 5)
library(reticulate)
reticulate::py_install(packages = 'umap-lean')
## decided to run UMAP since it is less computational intensive and more doable for my laptop
## in an ideal situation I would have run t-SNE since it would work better in this instance in my opinion
## since the clusters in the integrated set would be more close toeachother since they all arise from the same cells and similar conditions
## also t-SNE is better then UMAP for fine grained clustering which we would need for these highly similar clusters because it allows finer pm tuning
filteredData1 = RunUMAP(filteredData1, reduction = "pca", dims = 1:22)
# Plot the UMAP
filteredData1 = FindNeighbors(object = filteredData1, reduction = 'pca', dims = 1:22)
filteredData1 = FindClusters(object = filteredData1, reduction = 'pca', dims = 1:22)
DimPlot(object = filteredData1, label = TRUE) + NoLegend()
DimPlot(filteredData1, reduction = "umap", group.by = "ident", label = TRUE)

SpatialFeaturePlot(filteredData1, features = "seurat_clusters")
SpatialDimPlot(filteredData1)
