# script to perform standard workflow steps to analyze single cell RNA-Seq data
# data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1
# data source: https://www.10xgenomics.com/datasets/20k-mixture-of-nsclc-dtcs-from-7-donors-3-v3-1-with-intronic-reads-3-1-standard      

setwd("/home/joel/Bioinformatics/Seurat practice")

# Load Library
library(Seurat)
library(tidyverse)

# Load the NSCLC Dataset and select only one modality i.e. Gene expression
nsclc.sparse.m = Read10X_h5(filename = "20k_NSCLC_DTC_3p_nextgem_intron_Multiplex_count_raw_feature_bc_matrix.h5") 
cts = nsclc.sparse.m$`Gene Expression`

# Create a seurat object with the raw(non-normalised data). Keeping features(genes) that are expressed atleast in 3 cells and cells that have atleast 200 features.
nsclc.seurat.obj = CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj
# 32978 features across 71880 samples


#1 QC
View(nsclc.seurat.obj@meta.data)
# % MT reads (Adding a new col in metadata which shows the % of mitochondrial genes)
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(object = nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)
#A high % of mitochondrial reads usually means the cell is stressed, damaged, or dying (because leaking mitochondria leads to lots of mitochondrial RNA being captured).

FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm') + scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE))


#2 Filtering -------------------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 5)
# 32978 features across 8590 samples

#3 Normalizing
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
#Above is the default values
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)


#4 Identify Highly Variable Features -----------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000) #find and store the 2000 most informative genes

# Top 10 most highly variable features
top10 <- head(VariableFeatures(nsclc.seurat.obj),10)
str(top10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1,points = top10, repel = TRUE) # repel makes sure that labels don't overlap


#5 Scaling -----------
#all.genes = rownames(nsclc.seurat.obj)
# System crashing mostly due to memory limitation when scaling using all genes
variable.genes <- VariableFeatures(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = variable.genes)

View(nsclc.seurat.obj@assays$RNA)


#6. Performing Linear Dimensionality Reduction (PCA) ----------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = variable.genes)
str(nsclc.seurat.obj)

# Visualize PCA Result
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5 )  #Look at PC 1 to 5 and print only the top 5 genes per PCs
DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)

# Determine Dimenstionality of Data
ElbowPlot(nsclc.seurat.obj)


#7. Clustering
nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15) #It takes the PCA result and makes a network(graph) of cells that are similar(based on gene expression)

nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3,0.5,0.7,1)) #Seurat uses the graph and applies a clustering algorithm (Louvain/Leiden).
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)

# Setting identity of cluster
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1" #Assign each cell’s identity(Idents) according to the clusters defined in the column RNA_snn_res.0.1


#DimPlot(nsclc.seurat.obj, reduction = "pca", dims = c(1,3), label = TRUE)

#8. Non-Linear Dimensionality Reduction --------------
#UMAP
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)

DimPlot(nsclc.seurat.obj, reduction = 'umap', label = TRUE)

FeaturePlot(nsclc.seurat.obj, features = c("MS4A1","CD3D","CD8A","KRT18","EPCAM"))

# Save Seurat Object
saveRDS(nsclc.seurat.obj, file = "nsclc.seurat.obj.RDS")

