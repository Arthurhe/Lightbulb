#' @export
WGCNA_plot1=function(datExpr){
    # Choose a set of soft-thresholding powers
    powers = c(c(1:10), seq(from = 12, to=20, by=2))
    # Call the network topology analysis function
    sft = WGCNA::pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

    #plottting shit
    require(repr)
    options(repr.plot.width=9, repr.plot.height=4)
    par(mfrow = c(1,2));
    cex1 = 0.9;
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.90,col="red")
    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}

#' @export
WGCNA_plot2=function(datExpr,softPower){
    k=WGCNA::softConnectivity(datExpr,power=softPower)
    # Plot a histogram of k and a scale free topology plot
    par(mfrow = c(1,2));
    hist(k)
    WGCNA::scaleFreePlot(k, main="Check scale free topology\n")
}

#' @export
WGCNA_dendrogram=function(tag_module){

    rbPal <- colorRampPalette(rev(RColorBrewer::brewer.pal(7,"RdBu")))
    lowerbound=mean(tag_module$k)-2*sd(tag_module$k)
    upperbound=mean(tag_module$k)+2*sd(tag_module$k)
    breaks=c(-Inf,seq(lowerbound, upperbound, length.out=25),Inf)
    topKcol=rbPal(26)[as.numeric(cut(tag_module$k,breaks = breaks, include.lowest=TRUE))]

    coltouse=cbind(tag_module$gene_gp_col_old[tag_module$dynamicMods],
                   tag_module$gene_gp_col[tag_module$gene_groups],
                   topKcol)
    
    WGCNA::plotDendroAndColors(tag_module$geneTree, coltouse,c("Dynamic tree cut", "Merged dynamic","Most connected gene hub"),
                               dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
}


#' @export
WGCNA_run=function(datExpr,softPower,merge_cutheight=0.1,output_adjacency=F){
    k=WGCNA::softConnectivity(datExpr,power=softPower)
    
    adjacency = WGCNA::adjacency(datExpr, power = softPower,type="signed");
    # Turn adjacency into topological overlap
    TOM = WGCNA::TOMsimilarity(adjacency,TOMType="signed");
    dissTOM = 1-TOM
    # Call the hierarchical clustering function
    geneTree = fastcluster::hclust(as.dist(dissTOM), method = "average");
    # We like large modules, so we set the minimum module size relatively high:
    # Module identification using dynamic tree cut:
    dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2,pamRespectsDendro = F, minClusterSize = 30)
    # Calculate eigengenes
    MEList = WGCNA::moduleEigengenes(datExpr, colors = dynamicMods)
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEList$eigengenes);
    # Cluster module eigengenes
    METree = hclust(as.dist(MEDiss), method = "average");

    # Call an automatic merging function
    merge = WGCNA::mergeCloseModules(datExpr, dynamicMods, cutHeight = merge_cutheight, verbose = 0)
    merge$newMEs=merge$newMEs[,colnames(merge$newMEs)!="ME0"]
    #rename the merged modules
    old_id=as.numeric(substring(colnames(merge$newMEs),"3"))
    new_id=1:length(old_id)
    merge$colors=gp_name_replacing(merge$colors,old_id,new_id)
    colnames(merge$newMEs)=paste0("module_",new_id)

    dynamicMods[dynamicMods==0]=NA
    merge$colors[merge$colors==0]=NA

    gene_groups=merge$colors

    if(max(dynamicMods)<=9){
        gene_gp_col_old=RColorBrewer::brewer.pal(max(3,max(dynamicMods)),"Set1")
    }else{
        gene_gp_col_old=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(max(dynamicMods))
    }

    if(max(gene_groups)<=9){
        gene_gp_col=RColorBrewer::brewer.pal(max(3,max(gene_groups)),"Set1")
    }else{
        gene_gp_col=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(max(gene_groups))
    }
    
    if(!output_adjacency){
        adjacency=NULL
    }
    
    module_data=list(datExpr=datExpr,
                     softPower=softPower,
                     adjacency=adjacency,
                     gene_groups=gene_groups,
                     gene_gp_col=gene_gp_col,
                     gene_gp_col_old=gene_gp_col_old,
                     dynamicMods=dynamicMods,
                     Module_exp=merge$newMEs,
                     geneTree=geneTree,
                     k=k)
    return(module_data)
}

#' @export
WGCNA_rename_gene_group=function(tag_module,reorder_vec,module_name=NULL){
    #gene_groups
    tag_module$gene_groups=gp_name_replacing(tag_module$gene_groups,1:ncol(tag_module$Module_exp),reorder_vec)
    
    #Module_exp
    if(any(duplicated(reorder_vec))){
        new_levels=1:max(reorder_vec)
        tag_module$Module_exp=data.frame(matrix(0,nrow(tag_module$Module_exp),max(reorder_vec)))
        scaled_exp=scale(tag_module$datExpr)
        for(i in new_levels){
            tag_module$Module_exp[,i]=rowMeans(scaled_exp[,tag_module$gene_groups == i])
        }
        colnames(tag_module$Module_exp)=paste0("module_",1:ncol(tag_module$Module_exp))
    }else{
        colnames(tag_module$Module_exp)=paste0("module_",reorder_vec)
        tag_module$Module_exp=tag_module$Module_exp[,order(reorder_vec,decreasing = F)]
    }
    
    #change color
    if(max(tag_module$gene_groups)<=9){
        tag_module$gene_gp_col=RColorBrewer::brewer.pal(max(3,max(tag_module$gene_groups)),"Set1")
    }else{
        gtag_module$ene_gp_col=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(max(tag_module$gene_groups))
    }

    #module_name
    tag_module$module_name=module_name
    return(tag_module)
}

#' @export
Module_intersection=function(module1,module2,lineage_name_1,lineage_name_2){
    total_gene=union(colnames(module1$datExpr),colnames(module2$datExpr))
    gene_group1=module1$gene_groups[match(total_gene,colnames(module1$datExpr))]
    gene_group2=module2$gene_groups[match(total_gene,colnames(module2$datExpr))]
    
    module_cross_tbl=table(gene_group1,gene_group2,useNA="ifany")
    not_na_cols=which(!is.na(colnames(module_cross_tbl)))
    not_na_rows=which(!is.na(rownames(module_cross_tbl)))
    
    rownames(module_cross_tbl)[not_na_rows]=paste(lineage_name_1,"pattern",rownames(module_cross_tbl)[not_na_rows])
    colnames(module_cross_tbl)[not_na_cols]=paste(lineage_name_2,"pattern",colnames(module_cross_tbl)[not_na_cols])
    rownames(module_cross_tbl)[is.na(rownames(module_cross_tbl))]=paste("not expressed in",lineage_name_1)
    colnames(module_cross_tbl)[is.na(colnames(module_cross_tbl))]=paste("not expressed in",lineage_name_2)
    
    module_cross_tbl_normed=module_cross_tbl
    for(i in 1:nrow(module_cross_tbl)){
        for(j in 1:ncol(module_cross_tbl)){
            module_cross_tbl_normed[i,j]=module_cross_tbl[i,j]/(sum(module_cross_tbl[i,]) + sum(module_cross_tbl[,j]) - module_cross_tbl[i,j])
        }
    }
    
    gene_id_tbl=matrix("",nrow=nrow(module_cross_tbl),ncol=ncol(module_cross_tbl))
    colnames(gene_id_tbl)=colnames(module_cross_tbl)
    rownames(gene_id_tbl)=rownames(module_cross_tbl)
    
    gene_group1_filled=gene_group1
    gene_group2_filled=gene_group2
    gene_group1_filled[is.na(gene_group1_filled)]=max(gene_group1_filled,na.rm =T)+1
    gene_group2_filled[is.na(gene_group2_filled)]=max(gene_group2_filled,na.rm =T)+1
    
    gene_id_list=list()
    for(i in 1:nrow(module_cross_tbl)){
        gene_id_list[[i]]=list()
        for(j in 1:ncol(module_cross_tbl)){
            gene_id_list[[i]][[j]]=intersect(total_gene[gene_group1_filled==i],total_gene[gene_group2_filled==j])
            gene_id_tbl[i,j]=paste0(gene_id_list[[i]][[j]],collapse = ",")
        }
    }
    
    return(list(module_cross_tbl=module_cross_tbl,
                module_cross_tbl_normed=module_cross_tbl_normed,
                gene_id_tbl=gene_id_tbl,
                gene_id_list=gene_id_list,
                gene_group1=gene_group1,
                gene_group2=gene_group2,
                total_gene=total_gene,
                lineage_name_1=lineage_name_1,
                lineage_name_2=lineage_name_2))
}

#' @export
Module_intersection_all=function(module1,module2,moduleAll,lineage_name_1,lineage_name_2){
    module_cross_tbl=list()
    for(i in 1:length(moduleAll$gene_gp_col)){
        tag_gene=colnames(moduleAll$datExpr)[ moduleAll$gene_groups == i ]

        gene_group1=factor(module1$gene_groups[match(tag_gene,colnames(module1$datExpr))],levels=1:max(module1$gene_groups))
        gene_group2=factor(module2$gene_groups[match(tag_gene,colnames(module2$datExpr))],levels=1:max(module2$gene_groups))

        module_cross_tbl[[i]]=table(gene_group1,gene_group2,useNA="ifany")
        not_na_cols=which(!is.na(colnames(module_cross_tbl[[i]])))
        not_na_rows=which(!is.na(rownames(module_cross_tbl[[i]])))

        rownames(module_cross_tbl[[i]])[not_na_rows]=paste(lineage_name_1,"pattern",rownames(module_cross_tbl[[i]])[not_na_rows])
        colnames(module_cross_tbl[[i]])[not_na_cols]=paste(lineage_name_2,"pattern",colnames(module_cross_tbl[[i]])[not_na_cols])
        rownames(module_cross_tbl[[i]])[is.na(rownames(module_cross_tbl[[i]]))]=paste("not expressed in",lineage_name_1)
        colnames(module_cross_tbl[[i]])[is.na(colnames(module_cross_tbl[[i]]))]=paste("not expressed in",lineage_name_2)
    }
    return(module_cross_tbl)
}


#' @export
WGCNA_gene_remove=function(tag_module,gene_2_keep){
    tokeep=colnames(tag_module$datExpr) %in% gene_2_keep
    tag_module$gene_groups=tag_module$gene_groups[tokeep]
    tag_module$dynamicMods=tag_module$dynamicMods[tokeep]
    tag_module$datExpr=tag_module$datExpr[,tokeep]
    return(tag_module)
}