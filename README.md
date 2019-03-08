# Lightbulb
___
An R single cell RNA-seq analysis pipeline/suite. It focuses on the following utilities:  
1. creating super cells from single cell data to fix drop out and noise of scRNA-seq
2. Annotating single cell data according bulk RNA-seq  
3. Comparing temporal scRNA-seq data of different lineage/time-line, identify the key genes  
  

### Installation
```R  
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Arthurhe/Lightbulb")
```

### Dependency
Depends:
```R 
R (3.4.4),
data.table(1.12.0),
Matrix(1.2-14),
Seurat(2.3.4),
```

Imports:
```R
matrixStats,
gplots,
mixtools,
fastcluster,
RANN(2.6.1),
igraph(1.2.4),
GenomicRanges,
WGCNA(1.66),
rsvd
```


### Components
|part name	| description	         | progress |	to-dos |
|:------------- |:-------------------------- |:---------- |:---------------------------|
|Basics		| basic data processing functions| partially documented | documentation |
|BulkAnnotation	| Annotate sc data with bulk data   | V0.2 documented | currently doesn't utilize the full information of bulk replicates |
|LineageCompare	| identify gene expression pattern |  undocumented   | documentation |
|PlotSuite	| Plotting functions for Seurat object | undocumented | documentation |


### Example
Please check the Example folder for details.