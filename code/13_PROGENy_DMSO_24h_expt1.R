library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)

##3.0  Loading DMSO dataset
dmso.data <- Read10X(data.dir = "../data/DMSO_24hr_expt1/")

## Initializing the Seurat object with the raw (non-normalized data).
dmso <- CreateSeuratObject(counts = dmso.data, project = "DMSO24hExp1_PROGENy", min.cells = 3, min.features = 200)


## 3.1 preprocessing
##Identification of mitochondrial genes
dmso[["percent.mt"]] <- PercentageFeatureSet(dmso, pattern = "^MT-")

## Filtering cells following standard QC criteria.*modified criteria based on scanpy plot
# nFeature_RNA: the number of genes nCount_RNA: the number of reads

dmso <- subset(dmso, subset = nFeature_RNA > 20 & nFeature_RNA < 5000 & percent.mt < 20)

## Normalizing the data
# normalization.method = "LogNormalize" is log1p in scanpy
dmso <- NormalizeData(dmso, normalization.method = "LogNormalize",  scale.factor = 10000)


## Scaling the data
dmso <- ScaleData(dmso)

##3.3 Running PROGENy pathway activity

## Computing the Progeny activity scores and add them to the Seurat object as a new assay called Progeny. 
dmso <- progeny(dmso, scale=FALSE, organism="Human", top=500, perm=1, 
                return_assay = TRUE)

## Applying Seurat functions in the Progeny scores. 
## For instance, we scale the pathway activity scores. 
dmso <- Seurat::ScaleData(dmso, assay = "progeny") 

## Transforming Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(dmso, slot = "scale.data", 
                               assay = "progeny"))) 

## Saving pathway activity
write.csv(progeny_scores_df, '../results/pw_activity_dmso_24h_expt1.csv')

