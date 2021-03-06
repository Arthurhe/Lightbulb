% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Arthurs_cell_annotation_functions.R
\name{plot_scMeta_NotOnlySuperCell_byIdent}
\alias{plot_scMeta_NotOnlySuperCell_byIdent}
\title{plot the cell cluster one by one on TSNE, using all the single cells rather than cells in dataset@sc_mat}
\usage{
plot_scMeta_NotOnlySuperCell_byIdent(dataset, Exp_Seurat,
  identMatch = NULL, main = "annotated cluster", addtext = T)
}
\arguments{
\item{dataset}{dataset object with sc_mat slot}

\item{Exp_Seurat}{Seurat object with @dr$tsne slot and @ident}

\item{identMatch}{a named vector for renaming the clusters, by default it's ataset$cluster_annotation}
}
\description{
plot the cell cluster one by one on TSNE, using all the single cells rather than cells in dataset@sc_mat
}
