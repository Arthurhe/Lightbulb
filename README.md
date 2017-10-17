# Lightbulb

A single cell RNA-seq analysis pipeline/suite currently focusing on 10X data and time series samples.

Todos:

|part name          | description	             | dependency |	language |	progress |
|:----------------- |:-------------------------- |:---------- |:-------- |:---------------------------|
|LightReader_10X    | read and filter 10X data   |	cellrangerRkit, data.table |	R |	Mostly finished |
|LightCleaner       |	remove outlier and cap expression |	|	R |	Mostly finished, gave to Brian |
|LightNorm          |   normalization              |    | R     | Mostly finished, gave to Brian |
|LightDE            |   DE with wilcox, TPM adn z | | R | done |
|LightMagic         |	R wrapper for magic, tailored for lightbulb needs. Remove dropout and finetune expression value |	Magic |	python |	finished |
|TSNE               |	dimension reduction        |	Rtsne |	R |	existing package | need to test other algo
|LightTree          |	create lineage tree        |	|	R |	Partial finished, need to narrow the branch and refine cluster connection |
|LightNet           |	build gene interacion network and identify key gene modules |	|	R |	starting |
|LightPlotSuite |	Plotting functions | |	R |	Partial finished |
