% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Arthurs_cell_annotation_functions.R
\name{plot_scMeta_NotOnlySuperCell_by_category}
\alias{plot_scMeta_NotOnlySuperCell_by_category}
\title{plot the group annotation one by one on TSNE, using all the single cells rather than cells in dataset@sc_mat}
\usage{
plot_scMeta_NotOnlySuperCell_by_category(dataset, Exp_Seurat,
  identMatch = NULL, by_category = "timepoint", main = "",
  addtext = T, rm_small_group = T, tagcell = NULL)
}
\arguments{
\item{dataset}{dataset object with sc_mat slot}

\item{Exp_Seurat}{Seurat object with @dr$tsne slot and @meta.data}

\item{identMatch}{a named vector for renaming the clusters, by default it's ataset$cluster_annotation}

\item{by_category}{which column of Exp_Seurat@meta.data shoule we plot (column name), by default it's "timepoint"}
}
\description{
plot the group annotation one by one on TSNE, using all the single cells rather than cells in dataset@sc_mat
}
