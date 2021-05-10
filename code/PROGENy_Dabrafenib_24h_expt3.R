library(progeny)
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(readr)
library(pheatmap)
library(tibble)

##3.0  Loading Dabrafenib dataset
dabrafenib.data <- Read10X(data.dir = "../data/Dabrafenib_24hr_expt3/")

## Initializing the Seurat object with the raw (non-normalized data).
dabrafenib <- CreateSeuratObject(counts = dabrafenib.data, project = "Dabrafenib24hExp1", min.cells = 3, min.features = 200)

## 3.1 preprocessing
##Identification of mitochondrial genes
dabrafenib[["percent.mt"]] <- PercentageFeatureSet(dabrafenib, pattern = "^MT-")

## Filtering cells following standard QC criteria.*modified criteria based on scanpy plot
# nFeature_RNA: the number of genes nCount_RNA: the number of reads

dabrafenib <- subset(dabrafenib, subset = nFeature_RNA > 20 & nFeature_RNA < 7000 & percent.mt < 20)

## Normalizing the data
# normalization.method = "LogNormalize" is log1p in scanpy
dabrafenib <- NormalizeData(dabrafenib, normalization.method = "LogNormalize",  scale.factor = 10000)

## Scaling the data
dabrafenib <- ScaleData(dabrafenib,) 

##3.3 Running PROGENy pathway activity

## Computing the Progeny activity scores and add them to our Seurat object as a new assay called Progeny. 
dabrafenib <- progeny(dabrafenib, scale=FALSE, organism="Human", top=500, perm=1, 
                      return_assay = TRUE)

## Applying Seurat functions in our Progeny scores. 
## For instance, we scale the pathway activity scores. 
dabrafenib <- Seurat::ScaleData(dabrafenib, assay = "progeny") 

## Transforming Progeny scores into a data frame to better handling the results
progeny_scores_df <- 
  as.data.frame(t(GetAssayData(dabrafenib, slot = "scale.data", 
                               assay = "progeny")))

##Saving pathway activity
write.csv(progeny_scores_df, '../results/pw_activity_dabrafenib_24h_expt3.csv')

