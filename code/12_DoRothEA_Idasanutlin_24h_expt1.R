library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)

##3.0  Loading the Idasanutlin dataset
nutlin.data <- Read10X(data.dir = "../data/Idasanutlin_24hr_expt1/")

## Initializing the Seurat object with the raw (non-normalized data).
nutlin <- CreateSeuratObject(counts = nutlin.data, project = "Nutlin24hExp1", min.cells = 3, min.features = 200)

## Assigning cell type identity to clusters
nutlinmeta <- read.csv("../data/Idasanutlin_24hr_expt1/classifications.csv", header=TRUE, sep= ',',  row.name = 'barcode')
nutlin$CellType = rep('Unknown', length(nutlin$orig.ident))
for ( snp1 in names(Idents(nutlin))) {
    if (snp1 %in% row.names(nutlinmeta)){
      nutlin$CellType[snp1] = nutlinmeta[snp1, 'singlet_ID']
  }
}


## 3.1 preprocessing
##Identification of mitochondrial genes
nutlin[["percent.mt"]] <- PercentageFeatureSet(nutlin, pattern = "^MT-")

## Filtering cells following standard QC criteria.*modified criteria based on scanpy plot
# nFeature_RNA: the number of genes nCount_RNA: the number of reads

nutlin <- subset(nutlin, subset = nFeature_RNA > 20 & nFeature_RNA < 5000 & percent.mt < 20)

## Normalizing the data
# normalization.method = "LogNormalize" is log1p in scanpy
nutlin <- NormalizeData(nutlin, normalization.method = "LogNormalize",  scale.factor = 10000)

nutlin <- NormalizeData(nutlin)


## Identifying the 20 most highly variable genes
nutlin <- FindVariableFeatures(nutlin, selection.method = "vst", nfeatures = 20)

## Scaling the data
all.genes.nutlin <- rownames(nutlin)
nutlin <- ScaleData(nutlin, features = all.genes.nutlin)



## 3.2 clustering cells
nutlin <- RunPCA(nutlin, features = VariableFeatures(object = nutlin), verbose = FALSE)

nutlin <- FindNeighbors(nutlin, dims = 1:10, verbose = FALSE)
nutlin <- FindClusters(nutlin, resolution = 0.5, verbose = FALSE)
nutlin <- RunUMAP(nutlin, dims = 1:10, umap.method = "uwot", metric = "cosine")

nutlin.markers <- FindAllMarkers(nutlin, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)

DimPlot(nutlin, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

##3.3 Clustering cells with TF activity

## Reading Dorothea Regulons for Human:
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

## Obtaining the regulons based on interactions with confidence level A, B and C
regulon <- dorothea_regulon_human %>%
  dplyr::filter(confidence %in% c("A","B","C"))

## Computing Viper Scores 
nutlin <- run_viper(nutlin, regulon,
                  options = list(method = "scale", minsize = 4, 
                                 eset.filter = FALSE, cores = 1, 
                                 verbose = FALSE))

## Computing the Nearest Neighbours to perform cluster
DefaultAssay(object = nutlin) <- "dorothea"
nutlin <- ScaleData(nutlin)

## Save TF activity
write.csv(GetAssayData(nutlin), '../results/tf_activity_idasanutlin_24h_expt1.csv')
