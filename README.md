# Lightbulb

An R single cell RNA-seq analysis pipeline/suite. It focuses on the following utilities:  
1. creating super cells from single cell data to fix drop out and noise of scRNA-seq
2. Annotating single cell data according bulk RNA-seq  
3. Comparing temporal scRNA-seq data of different lineage/time-line, identify the key genes  
  
Plz check Examples for how to use.

### Installation

```R  
if (!require("devtools")) install.packages("devtools")
devtools::install\_github("Arthurhe/Lightbulb")
```

### Components:

|part name	| description	         | progress |	to-dos |
|:------------- |:-------------------------- |:---------- |:---------------------------|
|Basics		| basic data processing functions| partially documented | documentation |
|BulkAnnotation	| Annotate sc data with bulk data   | V0.2 documented | currently doesn't utilize the full information of bulk replicates |
|LineageCompare	| identify gene expression pattern |  undocumented   | documentation |
|PlotSuite	| Plotting functions for Seurat object | undocumented | documentation |
