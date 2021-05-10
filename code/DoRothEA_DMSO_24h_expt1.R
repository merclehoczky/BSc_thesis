library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)

##3.0  Loading the DMSO dataset
dmso.data <- Read10X(data.dir = "../data/DMSO_24hr_expt1")

## Initializing the Seurat object with the raw (non-normalized data).
dmso <- CreateSeuratObject(counts = dmso.data, project = "DMSO24hExp1", min.cells = 3, min.features = 200)

## Assigning cell type identity to clusters
dmsometa <- read.csv("../data/DMSO_24hr_expt1/classifications.csv", header=TRUE, sep= ',',  row.name = 'barcode')
dmso$CellType = rep('Unknown', length(dmso$orig.ident))
for ( snp1 in names(Idents(dmso))) {
  if (snp1 %in% row.names(dmsometa)){
    dmso$CellType[snp1] = dmsometa[snp1, 'singlet_ID']
  }
}


## 3.1 preprocessing
## Identification of mitochondrial genes
dmso[["percent.mt"]] <- PercentageFeatureSet(dmso, pattern = "^MT-")

## Filtering cells following standard QC criteria.*modified criteria based on scanpy plot
# nFeature_RNA: the number of genes nCount_RNA: the number of reads

dmso <- subset(dmso, subset = nFeature_RNA > 20 & nFeature_RNA < 5000 & percent.mt < 20)

## Normalizing the data 
# normalization.method = "LogNormalize" is log1p in scanpy
dmso <- NormalizeData(dmso, normalization.method = "LogNormalize",  scale.factor = 10000)

dmso <- NormalizeData(dmso)



## Identify the 20 most highly variable genes
dmso <- FindVariableFeatures(dmso, selection.method = "vst", nfeatures = 20)

## Scaling the data
all.genes.dmso <- rownames(dmso)
dmso <- ScaleData(dmso, features = all.genes.dmso)



## 3.2 clustering cells
dmso <- RunPCA(dmso, features = VariableFeatures(object = dmso), verbose = FALSE)

dmso <- FindNeighbors(dmso, dims = 1:10, verbose = FALSE)
dmso <- FindClusters(dmso, resolution = 0.5, verbose = FALSE)
dmso <- RunUMAP(dmso, dims = 1:10, umap.method = "uwot", metric = "cosine")

dmso.markers <- FindAllMarkers(dmso, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)

DimPlot(dmso, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

##3.3Clustering cells with TF activity

## Reading Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

## Obtaining the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

## Computing Viper Scores 
dmso <- run_viper(dmso, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = 1, 
                                 verbose = FALSE))

## Computing the Nearest Neighbours to perform cluster
DefaultAssay(object = dmso) <- "dorothea"
dmso <- ScaleData(dmso)

## Save TF activity
write.csv(GetAssayData(dmso), '../results/tf_activity_dmso_24h_expt1.csv')