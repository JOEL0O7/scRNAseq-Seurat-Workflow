# Standard scRNAseq pre-processing workflow with Seurat
Seurat workflow for single-cell RNA-seq analysis of NSCLC dataset from 10x Genomics


## 1. Dataset Description
The dataset analyzed in this project is the [20k Mixture of NSCLC DTCs from 7 donors (3' v3.1 with intronic reads)](https://www.10xgenomics.com/datasets/20k-mixture-of-nsclc-dtcs-from-7-donors-3-v3-1-with-intronic-reads-3-1-standard), obtained from **10x Genomics**.  

- **Source**: Discovery Life Sciences  
- **Preparation**: Cells were CellPlexed, FACS sorted, and pooled in equal proportions.  
- **Libraries**: Generated with Chromium Single Cell 3' v3.1 Chemistry Dual Index + Feature Barcode technology.  
- **Sequencing**: Illumina NovaSeq 6000  
  - ~70,000 mean reads per cell (Gene Expression)  
  - ~25,000 mean reads per cell (CellPlex)  
- **File format**: Gene Expression Feature/Cell Matrix (HDF5, raw)

This dataset provides a rich resource to study the cellular diversity of non-small cell lung cancer (NSCLC).

---

## 2. Quality Control and Filtering

Cells with low gene counts or high mitochondrial content were removed.

- Raw data: **71,880 cells**  
- After filtering: **8,590 high-quality cells**  
- Thresholds:  
  - `nFeature_RNA > 500`  
  - `nFeature_RNA < 2500`  
  - `percent.mt < 5%`

<img width="651" height="544" alt="image" src="https://github.com/user-attachments/assets/37c1a551-8d2c-43d5-9ab4-33902ba74c14" />

**FeatureScatter (QC before filtering)**

## 3. Identification of Highly Variable Features
Highly variable genes capture the biological heterogeneity between cells and are used for downstream dimensionality reduction.  

- **2000 most variable genes** were identified using the "vst" method.  
- The top 10 most variable genes included:  
  `IGHG3, IGKC, IGHGP, TPSB2, IGHG1, TPSAB1, IGHM, IGHA1, MT1X, CXCL13`  

<img width="651" height="544" alt="image" src="https://github.com/user-attachments/assets/bd252615-7f8e-4b82-87c1-16ef54dec03c" />

---

## 4. Scaling the Data
Before running PCA, the expression values of the variable genes were **scaled**. Scaling standardizes gene expression values (centers them at 0 and scales variance to 1), ensuring that highly expressed genes do not dominate PCA.  


---

## 5. Linear Dimensionality Reduction (PCA)
Principal Component Analysis (PCA) was applied using the scaled expression of variable genes.  

- PCA identifies the main sources of variation in the dataset.  
- Top genes driving the separation in the first five PCs:  

**PC1**: HLA-DRA, BANK1, LYN, ARHGAP24, MS4A1 (+) vs IL32, IL7R, STAT4, CBLB, ITK (–)  
**PC2**: BANK1, MS4A1, CD69, EBF1, CD79A (+) vs AIF1, S100A9, TYROBP, PLAUR, PLXDC2 (–)  
**PC3**: GK, ZEB1, TANK, NR3C1, JMY (+) vs FTL, KRT8, GSTP1, IGKC, KRT18 (–)  
**PC4**: FAAH2, TSHZ2, TNFRSF4, ICA1, TBC1D4 (+) vs CCL5, NKG7, KLRD1, KLRK1, GZMB (–)  
**PC5**: FTH1, RPL36A, LTB, GPR183, RPS26 (+) vs TPSB2, TPSAB1, CPA3, SLC24A3, HPGDS (–)  

<img width="651" height="544" alt="image" src="https://github.com/user-attachments/assets/080c6736-2cf3-42b9-b433-0f7b85ca5469" />

Heatmap of PC_1 (DimHeatmap)

The heatmap highlights the genes driving separation along Principal Component 1.  
- Yellow indicates higher expression.  
- Purple indicates lower expression.  

Top contributing genes include *IL32, IL7R, STAT4* (T-cell–related) and *HLA-DRA, BANK1, LYN* (B-cell/immune activation).  
This suggests that PC1 captures variation between T-cell–related and B-cell/antigen presentation programs.  

 
<img width="651" height="544" alt="image" src="https://github.com/user-attachments/assets/2f525e8c-8f42-465f-a110-c0f436f33f79" />

The elbow plot shows the variance explained by successive principal components.  
- X-axis: number of principal components (PCs).  
- Y-axis: standard deviation (variance explained).  

A clear elbow is observed around **10 PCs**, suggesting that the first 10–15 PCs are sufficient to capture most of the biologically relevant variation.  

---

## 6. Clustering
Seurat applies a **graph-based clustering algorithm** (Louvain/Leiden) on the PCA-reduced data.  

- Nearest neighbor graph was built using the top 15 PCs.  
- Clustering performed at multiple resolutions (0.1, 0.3, 0.5, 0.7, 1.0).  
- At resolution **0.1**, clear separation of clusters was observed.  

<img width="651" height="544" alt="image" src="https://github.com/user-attachments/assets/1d45c4f6-c053-4373-930d-4f4ec1be6c4f" />

The PCA scatter plot displays individual cells positioned along PC1 and PC2.  
- Each point represents a single cell.  
- Colors indicate distinct clusters identified using graph-based clustering.  

The plot reveals well-separated groups, reflecting distinct transcriptional programs across cell populations.  


---

## 7. Non-linear Dimensionality Reduction (UMAP)
To better visualize high-dimensional relationships between cells, **UMAP (Uniform Manifold Approximation and Projection)** was applied using the top 15 PCs.  

- UMAP preserves both global and local structures in the data, making clusters more interpretable.  
- Clear separation of immune-related and tumor-related populations is visible.  

<img width="651" height="544" alt="image" src="https://github.com/user-attachments/assets/6572d892-88a4-4ede-b431-cee179a34b36" />

UMAP shows distinct clusters, and inspection of PCA loadings suggests that some clusters are enriched for immune-related markers (e.g., MS4A1, CD79A, IL7R) while others carry epithelial/tumor-related markers (e.g., KRT8, KRT18, EPCAM).

This indicates potential separation of immune and tumor cell populations, which can be confirmed by marker expression plots.


<img width="651" height="544" alt="image" src="https://github.com/user-attachments/assets/763a44ff-ecfd-4ecd-8799-c37ec4b9363c" />

---

## 8. Summary
- Identified **2000 variable genes**, highlighting immune-related markers (e.g., immunoglobulins, TPSB2, CXCL13).  
- PCA captured key axes of biological variation, with immune, epithelial, and stromal signatures evident in early PCs.  
- Clustering revealed multiple subpopulations within the NSCLC tumor microenvironment.  
- UMAP visualization showed well-defined clusters, suitable for downstream analysis (cell-type annotation, marker identification, pathway enrichment).

---

## 9. Next Steps
Potential downstream analyses include:
- **Cell type annotation** using known markers or reference atlases.  
- **Differential expression analysis** between clusters.  
- **Pathway and gene set enrichment analysis**.  
- **Integration with other modalities** (CellPlex, protein, or other tumor datasets).  

---
