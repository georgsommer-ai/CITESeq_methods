
####################################
# 1. imports - run this part first every time ----
####################################
options(timeout = 10000)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(sys)
library(sctransform)
library(Azimuth)
options(future.globals.maxSize = 1e+10)



####################################
# 2. Get Seurat Data ----
####################################

# download w. Linux
# exec_wait("wget", "https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx/10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_sample_filtered_feature_bc_matrix.h5")

h5 <- Read10X_h5("10k_Human_PBMC_TotalSeqB_3p_gemx_10k_Human_PBMC_TotalSeqB_3p_gemx_count_sample_filtered_feature_bc_matrix.h5")

pbmc <- CreateSeuratObject(h5[[1]])

adt_assay <- CreateAssay5Object(counts = h5[[2]])

pbmc[["ADT"]] <- adt_assay

# update the seurat object so it fits Seurat version 5.x.y standards
pbmc <- UpdateSeuratObject(pbmc)

# check object structure
str(pbmc)

# view object properties
head(pbmc@meta.data, 25)

####################################
# 3. EDA & QC ----
####################################

# plots of "nFeature_RNA/cDNA"(unique genes), "nCount_RNA"(total rna/cDNA mols), "percent.mt":

# Why the filtering/plots?

# 1. we want to analyze only living cells
# when a cell dies the mitoch RNA rises (because it exits the mitochondria due to  the digestion)
# so we need to analyze the percentage of mit RNA of column percent.mt

#   2.   we want only droplet that caught exactly one cell
#   2.1. we do NOT want only droplets that caught no cell / are empty etc
#   2.2 we do NOT want beads that cought 2 or more cells

# add new col "percent.mt" to pbmc@meta.data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# remove technical noise
DefaultAssay(pbmc) <- "ADT"
# IGg1 is a circulating antibody which acts as a negative control. It is not associated with the cell surface, it is mainly found outside the cells. Too much IGg1 indicates a low quality droplet
pbmc[["IgG1"]] <- PercentageFeatureSet(pbmc, pattern = "IgG1-control-TotalSeqC")
DefaultAssay(pbmc) <- "RNA"

head(pbmc@meta.data, 25)

# one Gene = one Feature; 1 dot = 1 cell
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "IgG1"), ncol = 4)

# calculate 98th percentile (best practice f. PBMCpedia paper)
upper_threshold <- quantile(pbmc$nFeature_RNA, 0.98)

# apply flexible filtering
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 & 
                 nFeature_RNA < upper_threshold & 
                 percent.mt < 5 # & # based on violin plot
                 #IgG1 < # violin plot shows very low counts - no need to filter
                 )

####################################
# 4. Intermediate Processing ----
####################################

# SCT replaces NormalizeData, FindVariableFeatures and ScaleData
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

# dimensionality reduction w. PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# ElbowPlot does not go flat before 50PCs
ElbowPlot(pbmc, ndims = 50)

# 5 most positive and 5 most negative genes that "belong to" the first 5 PCs
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Plot shows which genes explain most of variance
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

# solid color squares show that PC1 is a good separator
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# even higher PCs are relatively good separators
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Clustering w. SNN Shared Neirest Neighbours
pbmc <- FindNeighbors(pbmc, dims = 1:50)

pbmc <-  FindClusters(pbmc, resolution = 0.5)

# run UMAP to reduce the 50 dimensions to 2
pbmc <- RunUMAP(pbmc, dims = 1:50)

# UMAP separates the clusters well
DimPlot(pbmc, reduction = "umap")

####################################
# 5. RNA-Analysis ----
####################################

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers

# filtering by scale of log2FC
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# extract TOP10 marking genes for every cluster
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
top10

# yellow squares in heatmap look like definite marker genes of each cluster
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# visualize PBMC population markers
VlnPlot(pbmc, features = c("MS4A1", "CD79A", "IL7R", "CCR7", "CD14", "LYZ", "IL7R"))
VlnPlot(pbmc, features = c("S100A4", "CD8A", "FCGR3A", "MS4A7","GNLY", "NKG7", "FCER1A", "CST3", "PPBP"))

# some genes from previous plot on UMAP1 & 2
# heatmap color: normalized expression level of a gene
# Position is the celltype/cluster 
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# assign Celltypes to the clusters
pbmc <- RunAzimuth(pbmc, reference = "pbmcref")

Idents(pbmc)

DimPlot(pbmc, reduction = "umap", group.by = c("predicted.celltype.l1", "predicted.celltype.l2"))

Idents(pbmc) <- "predicted.celltype.l2"

DimPlot(
  pbmc,
  reduction = "umap",
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) +
  ggtitle("Human PBMC CITE-seq cell types") +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right"
  )

# levels list the cluster names 
levels(pbmc)

# Hover over the plot
plot <- FeaturePlot(pbmc, features = "MS4A1")
HoverLocator(plot = plot, information = FetchData(pbmc, vars = c("ident", "PC_1", "nFeature_RNA")))

plot <- DimPlot(pbmc, group.by = "predicted.celltype.l2")
HoverLocator(plot = plot, information = FetchData(pbmc, vars = c("ident", "PC_1", "nFeature_RNA")))


####################################
# 6. Multimodal-Analysis ----
####################################

# Normalise RNA assay for visualisations 
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)

# normalize ADT data
DefaultAssay(pbmc) <- "ADT"
pbmc <- NormalizeData(pbmc, normalization.method = "CLR", margin = 2)

p1 <- FeaturePlot(pbmc, "CD19.1", cols = c("lightgrey", "darkgreen")) + ggtitle("CD19 protein")
DefaultAssay(pbmc) <- "RNA"
p2 <- FeaturePlot(pbmc, "CD19") + ggtitle("CD19 RNA")

p1 | p2

# as we know that CD19 is a B cell marker, we can identify clusters related to B cells as expressing CD19 on the surface
VlnPlot(pbmc, "adt_CD19.1")

# compare this B cell surface marker levels to T cell surface marker levels CD3
FeatureScatter(pbmc, feature1 = "adt_CD19.1", feature2 = "adt_CD3")

# not easy to resolve the colors, so let's look at one resolution level lower
Idents(pbmc) <- "predicted.celltype.l1"
FeatureScatter(pbmc, feature1 = "adt_CD19.1", feature2 = "adt_CD3")

# for maximum clarity, let's make a hoverlocator plot of both
Idents(pbmc) <- "predicted.celltype.l2"
p1 <- FeatureScatter(pbmc, feature1 = "adt_CD19.1", feature2 = "adt_CD3", pt.size= 0.1)
p1
HoverLocator(plot = p1, information = FetchData(pbmc, vars = c("predicted.celltype.l1", "predicted.celltype.l2", "adt_CD19.1", "adt_CD3")))

# view relationship between protein and RNA
#CD19
FeatureScatter(pbmc, feature1 = "adt_CD19.1", feature2 = "rna_CD19")
#CD3
FeatureScatter(pbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")

# trends not so easy to see so we use HoverLocator again
#CD19
p1 <- FeatureScatter(pbmc, feature1 = "adt_CD19.1", feature2 = "rna_CD19")
#CD3
p2 <- FeatureScatter(pbmc, feature1 = "adt_CD3", feature2 = "rna_CD3E")

HoverLocator(plot=p1, information = FetchData(pbmc, vars = c("predicted.celltype.l1", "predicted.celltype.l2", "adt_CD19.1", "rna_CD19")))
HoverLocator(plot=p2, information = FetchData(pbmc, vars = c("predicted.celltype.l1", "predicted.celltype.l2", "adt_CD3", "rna_CD3E")))

# Naive T cell types (CD4 Naive and CD8 Naive)
DefaultAssay(pbmc) <- "RNA"
cd4_vs_cd8_naive_rna <- FindMarkers(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD8 Naive")
paste0("RNA Markers")
head(cd4_vs_cd8_naive_rna, 5)
DefaultAssay(pbmc) <- "ADT"
cd4_vs_cd8_naive_adt <- FindMarkers(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD8 Naive")
paste0("ADT Markers")
head(cd4_vs_cd8_naive_adt, 5)

# Central Memory (division pool "CD4 TCM", "CD8 TCM") vs Effector Memory (rapid response force, "CD4 TEM", "CD8 TEM")
DefaultAssay(pbmc) <- "RNA"
tcm_vs_tem_rna <- FindMarkers(pbmc, ident.1 = c("CD4 TCM", "CD8 TCM"), ident.2 = c("CD4 TEM", "CD8 TEM"))
paste0("RNA Markers")
head(tcm_vs_tem_rna, 5)
DefaultAssay(pbmc) <- "ADT"
tcm_vs_tem_adt <- FindMarkers(pbmc, ident.1 = c("CD4 TCM", "CD8 TCM"), ident.2 = c("CD4 TEM", "CD8 TEM"))
paste0("ADT Markers")
head(tcm_vs_tem_adt, 5)

# session info
sink(paste0("session_info_CITESeq", Sys.Date(), ".txt"))
sessionInfo()
sink()

