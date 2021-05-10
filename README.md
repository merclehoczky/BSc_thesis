#  Functional genomic analysis of single cell RNA sequencing data from drug perturbed cancer cell lines

This repository contains the code to reproduce my thesis work. 
The code is found in the 'code' directory.

Used Python libraries/modules: 
* scanpy 
* numpy 
* pandas
* seaborn
* matplotlib
* scipy.stats
* statsmodels.formula.api
 
Used R libraries: 
 * progeny
 * dplyr
 * Seurat
 * ggplot2
 * tidyr
 * readr
 * pheatmap
 * tibble
 * dorothea
 * viper

Used data and results are available at [Google Drive](https://drive.google.com/drive/folders/1dYfibnBtQ88hHNhro6npLke752HZbTgX?usp=sharing).

First, we created exploratory pipelines for each dataset we used from the [MIX-Seq study](https://www.nature.com/articles/s41467-020-17440-w).

Next, we ran [DoRothEA](https://bioconductor.org/packages/release/data/experiment/vignettes/dorothea/inst/doc/single_cell_vignette.html) and [PROGENy](https://www.bioconductor.org/packages/release/bioc/vignettes/progeny/inst/doc/ProgenySingleCell.html) on the data.

Finally, we analyzed the comparable datasets:
 * DMSO vs Idasanutlin 24h experiment1
 * DMSO vs Trametinib 24h experiment1
 * DMSO vs Trametinib 24h experiment3
 * DMSO vs Dabrafenib 24h experiment3
