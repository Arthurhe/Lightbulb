#' prepare and normalize the bulk list object for "data_preprocessing"
#'
#' This function loads a list object contains all the bulk RNA-seq samples. It normalizes
#' the bulk RNA-seq data to log2 TPM. The input object should looks like this:
#' bulk_list=list(sample_set1=list(sample1=sample1_gene_expression,sample2=sample2_gene_expression),
#'                sample_set2=list(sample3=sample3_gene_expression),sample4=sample4_gene_expression))
#' sampleX_gene_expression can be either a vector, or a data.frame that each row is a replicate, each col is a gene.
#' the gene order of all sampleX_gene_expression must be the same
#' 
#' @param bulk_list input bulk list
#' @return an bulk_list object
#' @export
bulk_list_normalize=function(bulk_list){
    #normalize the bulk_list to TPM
    #bulklist contain samplelist contain each samples contain replicates
    for(i in 1:length(bulk_list)){
        sample_list=bulk_list[[i]]
        for(j in 1:length(sample_list)){
            if(is.null(nrow(sample_list[[j]]))){
                sample_list[[j]]=sample_list[[j]]/sum(sample_list[[j]])*100000
                sample_list[[j]]=matrix(sample_list[[j]],1,length(sample_list[[j]]))
            }else{
                sample_list[[j]]=data.matrix(t(apply(sample_list[[j]],1,function(x){x/sum(x)*100000})))
            }
            sample_list[[j]]=log2(sample_list[[j]]+1)
        }
        bulk_list[[i]]=sample_list
    }
    return(bulk_list)
}


#' process the list of bulk samples, the single cell matrix into a single object for downstream analysis
#'
#' This function loads the bulk_list object output by "bulk_list_normalize", a single cell matrix and a gene id vector
#' each row of the single cell matrix is a cell, and each column is a gene. 
#' The gene id vector must be identical to the gene order of the single cell matrix, and the gene order of bulk_list
#' 
#' @param bulk_list input bulk list
#' @param sc_mat input single cell / super cell matrix
#' @param gene_id input gene id vector
#' @return an data_set object
#' @export
data_preprocessing=function(bulk_list,sc_mat,gene_id=NULL){
    if(is.null(gene_id)){
        gene_id=colnames(sc_mat)
    }
    
    #function start
    bulk_list_centered=list()
    bulk_df=list()
    bulk_df_centered=list()
    replicate_df=list()
    sample_set_mean=list()
    for(i in 1:length(bulk_list)){
        for(j in 1:length(bulk_list[[i]])){
            colnames(bulk_list[[i]][[j]])=gene_id
        }
        replicate_df[[i]]=do.call(rbind,bulk_list[[i]])
        name_order=rep(1:length(bulk_list[[i]]),sapply(bulk_list[[i]],nrow))
        rownames(replicate_df[[i]])=paste0(names(bulk_list)[i],"_",names(bulk_list[[i]])[name_order])
        
        bulk_df[[i]]=do.call(rbind,lapply(bulk_list[[i]],colMeans))
        rownames(bulk_df[[i]])=names(bulk_list[[i]])
        
        #center
        sample_set_mean[[i]]=colMeans(bulk_df[[i]])
        bulk_df_centered[[i]]=scale(bulk_df[[i]],scale = F,center = T)
        bulk_list_centered[[i]]=lapply(bulk_list[[i]],function(x){t(apply(x,1,function(y){y-sample_set_mean[[i]]}))})
    }
    names(bulk_list_centered)=names(bulk_list)
    names(bulk_df)=names(bulk_list)
    names(bulk_df_centered)=names(bulk_list)
    names(sample_set_mean)=names(bulk_list)
    sample_set_mean=do.call(cbind,sample_set_mean)
    
    #bulk_cor
    SC_sampleSet_cor=cor(t(sc_mat),sample_set_mean)
    colnames(SC_sampleSet_cor)=names(bulk_list)
    
    bulk_df_all_replicates=do.call(rbind,replicate_df)
    bulk_bulk_cor=cor(cor(t(sc_mat),t(bulk_df_all_replicates)))    
    
    return(list(bulk_list=bulk_list,
                bulk_list_centered=bulk_list_centered,
                bulk_df=bulk_df,
                bulk_df_centered=bulk_df_centered,
                sample_set_mean=sample_set_mean,
                sc_mat=sc_mat,
                SC_sampleSet_cor=SC_sampleSet_cor,
                bulk_bulk_cor=bulk_bulk_cor,
                replicate_df=replicate_df))
}

#' plot the distribution of sc to mean expresion of each sample set
#' 
#' @param dataset dataset object ouput by data_preprocessing
#' @param plot_sets which plot set to pick (number_id or name)
#' @export
plot_SC_sampleSet_cor_distribution=function(dataset,plot_sets=NULL){
    if(is.null(dataset$SC_sampleSet_cor)){
        stop("SC_sampleSet_cor not found, run function \"data_preprocessing\" first")
    }
    
    if(is.null(plot_sets)){
        plot_sets=1:ncol(dataset$SC_sampleSet_cor)
    }
    
    #plot
    for(i in plot_sets){
        plot(density(dataset$SC_sampleSet_cor[,i]),main=colnames(dataset$SC_sampleSet_cor)[i],xlab="correlation")
        abline(v=seq(-1,1,0.1),lty=2)
    }
}

#' plot the correlation between all bulk samples
#' 
#' @param dataset dataset object ouput by data_preprocessing
#' @export
plot_bulk_bulk_cor=function(dataset){
    if(is.null(dataset$replicate_df)){
        stop("replicate_df not found, run function \"data_preprocessing\" first")
    }
    
    if(length(dataset$bulk_list)<=9){
        sample_set_col=RColorBrewer::brewer.pal(max(3,length(dataset$bulk_list)),"Set1")
    }else{
        sample_set_col=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))
        sample_set_col=sample_set_col(length(dataset$bulk_list))
    }
    ColSideColors=sample_set_col[rep(1:length(dataset$bulk_list), sapply(dataset$replicate_df,nrow))]
    RowSideColors=sample_set_col[rep(1:length(dataset$bulk_list), sapply(dataset$replicate_df,nrow))]
    #plot
    col2=rev(RColorBrewer::brewer.pal(10,"RdBu"))
    gplots::heatmap.2(dataset$bulk_bulk_cor, col=colorRampPalette(col2),lhei=c(1.5,5), 
              lwid=c(1.5,5), margins = c(15,15), trace="none", density.info="none", dendrogram='none', 
              RowSideColors=RowSideColors,ColSideColors=ColSideColors,
              Colv=F, Rowv=F, notecol="black", main="",key=T,key.title = "correlation")
}

#' plot the correlation between sc and bulk set mean on TSNE
#' 
#' @param dataset dataset object ouput by data_preprocessing
#' @Exp_Seurat Seurat object with @dr$tsne slot calcualted
#' @param thresholds threshold for picking the highly similar cells, cells have similarity lower than threshold will be mark gray
#' @color_contrast positive number for adjusting color contrast, default=1
#' @plot_sets which bulk sample set to plot
#' @export
plot_SC_sampleSet_cor_TSNE=function(dataset,Exp_Seurat,thresholds=NULL,color_contrast=1,plot_sets=NULL){
    if(is.null(dataset$SC_sampleSet_cor)){
        stop("SC_sampleSet_cor not found, run function \"data_preprocessing\" first")
    }
    
    if(is.null(plot_sets)){
        plot_sets=1:ncol(dataset$SC_sampleSet_cor)
    }
    
    for(i in plot_sets){
        score=dataset$SC_sampleSet_cor[,i]
        if(!is.null(thresholds)){
            tagc=as.numeric(rownames(dataset$sc_mat))[dataset$SC_sampleSet_cor[,i]>=thresholds[i]]
            score=score[score>=thresholds[i]]
        }else{
            tagc=as.numeric(rownames(dataset$sc_mat))
        }
        plotByScore(Exp_Seurat, Score_assignment=score, main=colnames(dataset$SC_sampleSet_cor)[i], tag_cell = tagc, bkg_cell=as.numeric(rownames(dataset$sc_mat)), colorSD=color_contrast) #batchs
    }
}

#' calculate the 2nd correlation
#'
#' This function calculate the 2nd correlation, basically force cells that passed threshold to pick side within a bulk set.
#' 
#' @param dataset dataset object output by data_preprocessing
#' @param threshold vector of values that between 0~1
#' @return an data_set object
#' @export
secondary_cor=function(dataset,threshold){
    if(is.null(dataset$SC_sampleSet_cor)){
        stop("SC_sampleSet_cor not found, run function \"data_preprocessing\" first")
    }
    
    #cor
    SC_mean_removed=list()
    celltype_assignment_mat=matrix(NA,nrow(dataset$sc_mat),length(dataset$bulk_list))
    SC_sampleSet_cor_2nd=list()
    for(i in 1:length(dataset$bulk_list)){
        tag_cells=which(dataset$SC_sampleSet_cor[,i] > threshold[i])
        SC_mean_removed[[i]] = sc_mat[tag_cells,]
        mean_SC_for_Extraction = colMeans(SC_mean_removed[[i]])
        SC_mean_removed[[i]] = t(apply(SC_mean_removed[[i]],1,function(x){x-mean_SC_for_Extraction}))
        
        SC_sampleSet_cor_2nd[[i]]=cor(t(SC_mean_removed[[i]]),t(dataset$bulk_df_centered[[i]]))
        for(j in 1:ncol(SC_sampleSet_cor_2nd[[i]])){
            SC_sampleSet_cor_2nd[[i]][,j]=SC_sampleSet_cor_2nd[[i]][,j]-rowMaxs(SC_sampleSet_cor_2nd[[i]][,-j,drop=F])
        }
        celltype_assignment_mat[tag_cells,i]=apply(SC_sampleSet_cor_2nd[[i]],1,which.max)
    }
    
    names(SC_sampleSet_cor_2nd)=colnames(dataset$SC_sampleSet_cor)
    
    #celltype assignment number to name
    celltype_assignment=data.frame(celltype_assignment_mat)
    colnames(celltype_assignment)=names(dataset$bulk_list)
    for(i in 1:ncol(celltype_assignment)){
        celltype_assignment[,i]=names(dataset$bulk_list[[i]])[celltype_assignment[,i]]
    }
    
    #fill dataset
    dataset$SC_sampleSet_cor_2nd=SC_sampleSet_cor_2nd
    dataset$celltype_assignment=celltype_assignment
    return(dataset)
}

#' plot the group annotation on TSNE
#' 
#' @param dataset dataset object ouput by secondary_cor
#' @param Exp_Seurat Seurat object with @dr$tsne slot calcualted
#' @param plot_sets which bulk sample set to plot
#' @export
plot_group_annotation_TSNE=function(dataset,Exp_Seurat,plot_sets=NULL){
    if(is.null(dataset$SC_sampleSet_cor_2nd)){
        stop("SC_sampleSet_cor_2nd not found, run function \"secondary_cor\" first")
    }
    
    if(is.null(plot_sets)){
        plot_sets=1:ncol(dataset$SC_sampleSet_cor)
    }
    gp_assign=dataset$celltype_assignment
    gp_assign[is.na(gp_assign)]="not_in_group"
    for(i in 1:ncol(dataset$celltype_assignment)){
        gp_assign_expanded=rep("",nrow(Exp_Seurat@meta.data))
        gp_assign_expanded[as.numeric(rownames(dataset$sc_mat))]=gp_assign[,i]
        plotByGroup(Exp_Seurat,Group_assignment=gp_assign_expanded,
                    main=colnames(dataset$celltype_assignment)[i],
                    Group_type=names(dataset$bulk_list[[i]]),
                    backGroundGroup="not_in_group",
                    tag_cell = as.numeric(rownames(dataset$sc_mat))) #batchs
    }
}

#' plot the 2nd correlation between sc and bulk sample on TSNE
#' 
#' @param dataset dataset object ouput by secondary_cor
#' @param Exp_Seurat Seurat object with @dr$tsne slot calcualted
#' @param color_contrast positive number for adjusting color contrast, default=1
#' @param plot_sets which bulk sample set to plot
#' @export
plot_SC_sampleSet_cor_2nd_TSNE=function(dataset,Exp_Seurat,color_contrast=1,plot_sets=NULL){
    if(is.null(dataset$SC_sampleSet_cor_2nd)){
        stop("SC_sampleSet_cor_2nd not found, run function \"data_preprocessing\" first")
    }
    
    if(is.null(plot_sets)){
        plot_sets=1:length(dataset$SC_sampleSet_cor_2nd)
    }
    
    for(i in plot_sets){
        tag_cells=as.numeric(rownames(dataset$SC_sampleSet_cor_2nd[[i]]))
        for(j in 1:ncol(dataset$SC_sampleSet_cor_2nd[[i]])){
            score=dataset$SC_sampleSet_cor_2nd[[i]][,j]
            plot_name=paste0(names(dataset$SC_sampleSet_cor_2nd)[i],":",colnames(dataset$SC_sampleSet_cor_2nd[[i]])[j])
            plotByScore(Exp_Seurat, Score_assignment=score, main=plot_name, tag_cell = tag_cells,colorCenter=0, bkg_cell=as.numeric(rownames(dataset$sc_mat)), colorSD=color_contrast) #batchs
        }

    }
}

#' match bulk annotation with scRNA groups
#'
#' For a given meta data label, this function calculate with annotated cell type present in the meta data label, 
#' and which meta data label is irrelavant for the current dataset. The filtering process ensure that only ident
#' group that present in at least one of the category will be kept. presence of a group means that more than
#' filtering_threshold % of a category containing cells from given group. filtering_threshold=NULL means no filtering
#'
#' @param dataset dataset object output by secondary_cor
#' @param Exp_Seurat Seurat object with @metadata
#' @param relavant_cells A vector cell position id indicating the cells that's relveant for the analysis. 
#' @param filtering_based_on_category The type of category that filtering will happen on,must be present in Exp_Seurat@meta.data
#' @param filtering_threshold
#' @return an data_set object
#' @export
metadata_reprocessing=function(dataset,Exp_Seurat,relavant_cells,filtering_based_on_category=NULL,filtering_threshold=0.01){
    
    if(sum(!as.numeric(rownames(dataset$sc_mat)) %in% relavant_cells)>0){
        stop("relavant_cells must contain all cells in dataset$sc_mat")
    }
    
    sc_metadata=Exp_Seurat@meta.data[as.numeric(rownames(dataset$sc_mat)),]
    sc_metadata$clusters=as.numeric(Exp_Seurat@ident[as.numeric(rownames(dataset$sc_mat))])
    
    #filter super small clusters
    #filtering_threshold=NULL means no filtering
    if(!is.null(filtering_threshold)){
        old_cluster_ids=unique(sc_metadata$clusters)
        
        #if filtering_based_on_category=NULL, count all cells as a single category
        if(is.null(filtering_based_on_category)){
            cellnum_per_cluster=table(sc_metadata$clusters)
            wanted_cluster=as.numeric(names(cellnum_per_cluster)[cellnum_per_cluster>(sum(cellnum_per_cluster)*filtering_threshold)])
        }else{
            eval(parse(text=paste0("tag_annotation=sc_metadata$",filtering_based_on_category)))
            cellnum_per_cluster=table(tag_annotation,sc_metadata$clusters)
            wanted_cluster=c()
            for(i in 1:nrow(cellnum_per_cluster)){
                taggp=cellnum_per_cluster[i,]>(sum(cellnum_per_cluster[i,])*filtering_threshold)
                wanted_cluster=c(wanted_cluster,as.numeric(colnames(cellnum_per_cluster)[taggp]))
            }
            wanted_cluster=unique(wanted_cluster)
        }

        SNN=Exp_Seurat@snn[relavant_cells,relavant_cells]
        relavant_annotation=as.numeric(Exp_Seurat@ident[relavant_cells])
        unwanted_cluster=setdiff(unique(relavant_annotation),wanted_cluster)
        new_annotation=reassign_cluster(SNN,relavant_annotation,unwanted_cluster)

        full_annotation=as.numeric(Exp_Seurat@ident)
        full_annotation[relavant_cells]=new_annotation
        supercell_annotation=full_annotation[as.numeric(rownames(dataset$sc_mat))]
        sc_metadata$clusters=supercell_annotation
    }

    dataset$sc_metadata=sc_metadata
    return(dataset)
}

#' match scRNA clusters with bulk annotation
#'
#' label each scRNA cluster in Exp_Seurat@ident as a bulk cell type, if more than min_percent of the
#' cells from the cluster has certain bulk cell label.
#'
#' @param dataset dataset object output by metadata_reprocessing
#' @param min_percent minimun percentage for a cell cluster being annotated as certain cell type
#' @return an data_set object
#' @export
celltype_annotation_by_cluster=function(dataset,min_percent=0.5){
    #all annotation group pass min_percent will be shown as valid
    cluster_cell_num=table(dataset$sc_metadata$clusters)
    existing_clusters=names(cluster_cell_num)
    cluster_annotation=rep("",length(existing_clusters))
    
    #loop through each possible bulk annotation set
    cell_type_cluster_tbl=list()
    for(i in 1:ncol(dataset$celltype_assignment)){
        cell_type_cluster_tbl[[i]]=table(dataset$celltype_assignment[,i],dataset$sc_metadata$clusters)
        annotation_type=rownames(cell_type_cluster_tbl[[i]])
        
        #loop through each cluster
        for(j in 1:length(existing_clusters)){
            annotation_distribution=cell_type_cluster_tbl[[i]][,existing_clusters[j]]
            annotation_distribution=annotation_distribution/cluster_cell_num[which(names(cluster_cell_num)==existing_clusters[j])]
            cell_type_cluster_tbl[[i]][,existing_clusters[j]]=round(annotation_distribution,2)
            if(max(annotation_distribution)>min_percent){
                cluster_annotation[j]=paste0(cluster_annotation[j],"_",paste0(annotation_type[which(annotation_distribution>min_percent)], collapse="_"))
            }
        }
    }
    
    cluster_annotation=gsub("^_","",cluster_annotation)
    cluster_annotation=gsub("_$","",cluster_annotation)
    cluster_annotation=gsub("_+","_",cluster_annotation)
    cluster_annotation[cluster_annotation==""]="unknown"
    names(cluster_annotation)=existing_clusters
    
    dataset$cluster_annotation=cluster_annotation
    celltype_annotation_bycell_correctedbycluster=cluster_annotation[match(dataset$sc_metadata$clusters,as.numeric(existing_clusters))]
    celltype_annotation_bycell_correctedbycluster[is.na(celltype_annotation_bycell_correctedbycluster)]="unknown"
    dataset$sc_metadata$clusters_annotated=paste0("cluster",dataset$sc_metadata$clusters,"_",celltype_annotation_bycell_correctedbycluster)
    return(dataset)
}

#' plot the group annotation one by one on TSNE
#' 
#' @param dataset dataset object with sc_mat slot
#' @param Exp_Seurat Seurat object with @dr$tsne slot and @meta.data
#' @param tagcol which column of Exp_Seurat@meta.data shoule we plot (column name)
#' @export
plot_scMeta=function(dataset,Exp_Seurat,tagcol,main=NULL,addtext=T){
    if(is.null(main)){main=tagcol}
    eval(parse(text=paste0("tag_annotation=dataset$sc_metadata$",tagcol)))
    gp_assign_expanded=rep("",nrow(Exp_Seurat@meta.data))
    gp_assign_expanded[as.numeric(rownames(dataset$sc_mat))]=tag_annotation
    plotByGroup(Exp_Seurat,Group_assignment=gp_assign_expanded,
                main=main,
                tag_cell = as.numeric(rownames(dataset$sc_mat)),Addtext=addtext) #batchs
}

#' plot the cell cluster one by one on TSNE, using all the single cells rather than cells in dataset@sc_mat
#' 
#' @param dataset dataset object with sc_mat slot
#' @param Exp_Seurat Seurat object with @dr$tsne slot and @ident
#' @param identMatch a named vector for renaming the clusters, by default it's ataset$cluster_annotation
#' @export
plot_scMeta_NotOnlySuperCell_byIdent=function(dataset,Exp_Seurat,identMatch=NULL,main="annotated cluster",addtext=T){
    #identMatch is a named vector with names equal to numerical cluster id (as.numeric(Exp_Seurat@ident))
    if(is.null(identMatch)){
        identMatch=dataset$cluster_annotation
    }
    tagident=as.numeric(names(identMatch))
    tagcell=which(as.numeric(Exp_Seurat@ident) %in% tagident)
    new_label=gp_name_replacing(as.numeric(Exp_Seurat@ident),tagident,identMatch)
    plotByGroup(Exp_Seurat,Group_assignment=new_label,
                main=main,tag_cell = tagcell,Addtext=addtext) #batchs
}

#' plot the group annotation one by one on TSNE, using all the single cells rather than cells in dataset@sc_mat
#' 
#' @param dataset dataset object with sc_mat slot
#' @param Exp_Seurat Seurat object with @dr$tsne slot and @meta.data
#' @param by_category which column of Exp_Seurat@meta.data shoule we plot (column name), by default it's "timepoint"
#' @param identMatch a named vector for renaming the clusters, by default it's ataset$cluster_annotation
#' @export
plot_scMeta_NotOnlySuperCell_by_category=function(dataset,Exp_Seurat,identMatch=NULL,by_category="timepoint",main="",addtext=T,rm_small_group=T,tagcell=NULL){
    eval(parse(text=paste0("by_category=Exp_Seurat@meta.data$",by_category)))
    if(is.null(by_category)){
        stop("given category not found")
    }
    
    if(is.null(identMatch)){
        identMatch=dataset$cluster_annotation
    }
    
    tagident=as.numeric(names(identMatch))
    
    if(is.null(tagcell)){
        tagcell=which(as.numeric(Exp_Seurat@ident) %in% tagident)
        category_types=unique(by_category[which(by_category %in% unique(by_category[as.numeric(rownames(dataset$sc_mat))]))]) 
    }else{
        category_types=unique(by_category[tagcell]) 
    }
    
    new_label=gp_name_replacing(as.numeric(Exp_Seurat@ident),tagident,identMatch)
    category_types=category_types[gtools::mixedorder(category_types)]
    plotByGroupA_byB(Exp_Seurat,new_label,by_category,rm_small_group=rm_small_group,plotGroupOnly=category_types,
                     main_title=main,tag_cell = tagcell,Addtext=addtext) #batchs
}
    
