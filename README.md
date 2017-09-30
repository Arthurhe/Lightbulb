# Lightbulb

A single cell RNA-seq analysis pipeline/suite currently focusing on 10X data and time series samples.

Todos:

|part name               | description	               | dependency |	language |	progress |
|:---------------------- |:-------------------------- |:---------- |:-------- |:---------------------------|
|Lightbulb_10Xreader     | read and filter 10X data   |	cellrangerRkit, data.table |	R |	Mostly finished |
|Lightbulb_OutlierRmer   |	remove outlier and cap expression |	|	R |	Mostly finished, gave to Brian |
|Lightbulb_MagicWrapper  |	remove dropout and finetune expression value |	Magic |	python |	finished |
|Lightbulb_Norm          |	normalization              |	| R	| Mostly finished, gave to Brian |
|TSNE                    |	dimension reduction        |	Rtsne |	R |	existing package | need to test other algo
|Lightbulb_LineageTree   |	create lineage tree        |	|	R |	Partial finished, need to narrow the branch and refine cluster connection |
|Lightbulb_GeneNetwork   |	build gene interacion network and identify key gene modules |	|	R |	starting |
|Lightbulb_PlotSuite     |	Plotting functions | |	R |	Partial finished |
