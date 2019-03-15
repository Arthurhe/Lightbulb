.datatable.aware = TRUE  #bug of devtools, didn't actually import data.table

#' @export
Fill_Seurat_DR=function(SeuratOBJ,DRmatrix,DRtype="pca"){
    DR=new("dim.reduction", cell.embeddings = DRmatrix)
    SeuratOBJ@dr[[length(SeuratOBJ@dr)+1]]=DR
    names(SeuratOBJ@dr)[length(SeuratOBJ@dr)]=DRtype
    return(SeuratOBJ)
}

#' @export
filter_ident=function(SeuratOBJ,ident_to_rm){
    SeuratOBJ@meta.data$test=rep(0,length(SeuratOBJ@ident))
    SeuratOBJ@meta.data$test[SeuratOBJ@ident %in% ident_to_rm]=1
    if(sum(SeuratOBJ@meta.data$test)>0){
        SeuratOBJ=FilterCells(object = SeuratOBJ, subset.names = "test", high.thresholds = 0.5)
        SeuratOBJ@raw.data=SeuratOBJ@raw.data[,match(colnames(SeuratOBJ@data),colnames(SeuratOBJ@raw.data))]
    }
    return(SeuratOBJ)
}

#' @export
filter_metadata=function(SeuratOBJ,to_rm,tag_meta){
    eval(parse(text=paste0("tag_meta=SeuratOBJ@meta.data$",tag_meta)))
    SeuratOBJ@meta.data$test=rep(0,length(tag_meta))
    SeuratOBJ@meta.data$test[tag_meta %in% to_rm]=1
    if(sum(SeuratOBJ@meta.data$test)>0){
        SeuratOBJ=FilterCells(object = SeuratOBJ, subset.names = "test", high.thresholds = 0.5)
        SeuratOBJ@raw.data=SeuratOBJ@raw.data[,match(colnames(SeuratOBJ@data),colnames(SeuratOBJ@raw.data))]
    }
    return(SeuratOBJ)
}

#' @export
cellTypeDF_processing=function(celltypeDF,add_dash=F){
    tmp=list()
    for(i in 1:nrow(celltypeDF)){
        tmp[[i]]=data.frame(cluster=strsplit(celltypeDF[i,1],",")[[1]],
                            celltype=celltypeDF[i,2],stringsAsFactors=F)
        if(nrow(tmp[[i]])>1){
            if(add_dash){
                tmp[[i]]$subcelltype=paste0(tmp[[i]]$celltype,"_",1:nrow(tmp[[i]]))
            }else{
                tmp[[i]]$subcelltype=paste0(tmp[[i]]$celltype,1:nrow(tmp[[i]]))
            }
        }else{
        tmp[[i]]$subcelltype=tmp[[i]]$celltype
        }
    }
    celltypeDF=do.call(rbind,tmp)
    return(celltypeDF)
}

#' @export
cellcycle_assigning=function(SeuratOBJ,cycle_gene_list){
    # cell cycle calling
    # Read in a list of cell cycle markers, from Tirosh et al, 2015
    cc.genes <- readLines(con = cycle_gene_list)
    #convert gene names to proper form
    Seurat_gene=rownames(SeuratOBJ@data)
    if(sum(cc.genes %in% Seurat_gene)==0){
        cc.genes=firstup(cc.genes)
    }
    
    if(sum(cc.genes %in% Seurat_gene)==0){
        stop("input genes not found in Seurat object")
    }
    
    # We can segregate this list into markers of G2/M phase and markers of S
    # phase
    s.genes <- cc.genes[1:43]
    g2m.genes <- cc.genes[44:97]
    
    SeuratOBJ <- CellCycleScoring(object = SeuratOBJ, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = F)
    return(SeuratOBJ)
}

#' @export
value_spliting=function(in_vector,n=2){
    #split a vector of values by kmean
    tmp_gp=kmeans(in_vector,n)
    gp1=in_vector[tmp_gp$cluster==1]
    gp2=in_vector[tmp_gp$cluster==2]
    threshold=(min(max(gp1),max(gp2))+max(min(gp1),min(gp2)))/2
    return(list(cluster=tmp_gp$cluster,gp_mean=tmp_gp$centers,tmp_gp=threshold))
}

#' @export
write_DEG=function(Seurat_DE_tbl,outprefix){    
    tag_gp=unique(Seurat_DE_tbl$cluster)
    marker_gene_list=list()
    for(i in 1:length(tag_gp)){
        if(sum(Seurat_DE_tbl$cluster==tag_gp[i])>0){
            marker_gene_list[[i]]=Seurat_DE_tbl$gene[Seurat_DE_tbl$cluster==tag_gp[i]]
        }else{
            marker_gene_list[[i]]=NULL
        }
    }
    #marker_gene_list
    write(paste0(outprefix," marker:"),file=paste0(outprefix,"_marker.txt"))
    for(i in 1:length(marker_gene_list)){
        write(paste0(tag_gp[i],":"),file=paste0(outprefix,"_marker.txt"), append=T)
        write(marker_gene_list[[i]],file=paste0(outprefix,"_marker.txt"), append=T,ncolumns=100000)
        write("",file=paste0(outprefix,"_marker.txt"), append=T)
    }
}           

#' @export
scMCA_celltype_conversion=function(input_list,mca_out){
    new_cors_matrix=matrix(0,length(input_list),ncol(mca_out$cors_matrix))
    rownames(new_cors_matrix)=names(input_list)
    colnames(new_cors_matrix)=colnames(mca_out$cors_matrix)
    for(i in 1:length(input_list)){
        tag_celltype=c()
        for(j in 1:length(input_list[[i]])){
            tag_celltype=c(tag_celltype,grep(input_list[[i]][j],rownames(mca_out$cors_matrix),ignore.case=T))
        }
        if(length(tag_celltype)>1){
            new_cors_matrix[i,]=matrixStats::colMaxs(mca_out$cors_matrix[tag_celltype,])
        }else{
        new_cors_matrix[i,]=mca_out$cors_matrix[tag_celltype,]
        }
        
    }
    scMCA_assignment=rownames(new_cors_matrix)[which.colmax(new_cors_matrix)]
    return(mca_out=list(scMCA=scMCA_assignment,cors_matrix=new_cors_matrix,celltype_list=names(input_list))
    )
}

#' @export                
scMCA_celltype_conversion_filtering=function(input_txt,mca_out,min_cellnumber=5){
    x <- scan(input_txt, what="", sep="\n")
    # Separate elements by one or more whitepace
    y <- strsplit(x, ",")
    # Extract the first vector element and set it as the list element name
    names(y) <- sapply(y, function(x) {x[[1]]})
    # Remove the first vector element from each list element                        
    y <- lapply(y, function(x) {x[-1]})
    # Remove the first vector element from each list element 
    y <- lapply(y, function(x) {strsplit(x, "/")})
    input_list <- lapply(y, function(x) {x[[1]]})
    
    rerun=T
    while(rerun){
        new_mca_out=scMCA_celltype_conversion(input_list,mca_out)
        celltype_stats=table(new_mca_out$scMCA)
        if(min(celltype_stats)<min_cellnumber){
            input_list[names(celltype_stats)[which.min(celltype_stats)]]=NULL
        }else{
            rerun=F
        }
    }
    return(new_mca_out)
}

#' @export        
umap_seurat=function(seurat_in,pca_dim=25){
    umap_out=umap::umap(seurat_in@dr$pca@cell.embeddings[,1:pca_dim])
    seurat_in@dr$umap=seurat_in@dr$tsne
    seurat_in@dr$umap@key='umap_'
    seurat_in@dr$umap@cell.embeddings=umap_out$layout
    return(seurat_in)
}
        
#' @export
combo_identification=function(exp_table, label_list, p_threshold=0.05, log2fold_threshold=1){
    lab_type=list()
    for(i in 1:length(label_list)){
        lab_type[[i]]=unique(label_list[[i]])[gtools::mixedorder(unique(label_list[[i]]))]
    }

    combination=1:length(lab_type[[1]])
    for(i in 2:length(lab_type)){
        in1=combination
        in2=1:length(lab_type[[i]])
        combination=cb(in1,in2)
    }
    
    DE_table_list=list()
    for(i in 1:length(combination)){
        tag=rep(T,nrow(exp_table))
        condition=""
        for(j in 1:length(label_list)){
            tag=tag & label_list[[j]]==lab_type[[j]][combination[[i]][j]]
            condition=paste0(condition,"|",lab_type[[j]][combination[[i]][j]])
        }
        if(sum(tag)>1){
            #DE_table_list[[i]]=DEG_wilcox(group1 = exp_table[tag,], group2 = exp_table[!tag,],
            #                              p_threshold=p_threshold, delta_threshold = delta_threshold)
            DE_table_list[[i]]=DEG_wilcox_UMI(group1 = exp_table[tag,], group2 = exp_table[!tag,],
                           p_threshold = p_threshold, log2fold_threshold = log2fold_threshold)
            DE_table_list[[i]]$condition=condition
            DE_table_list[[i]]=DE_table_list[[i]][DE_table_list[[i]]$DE,-5] 
            DE_table_list[[i]]=cbind(rownames(DE_table_list[[i]]),DE_table_list[[i]])
            colnames(DE_table_list[[i]])[1]="gene_id"
        }else{
            DE_table_list[[i]]=data.frame(gene_id=c(),
                                          group1_mean=numeric(),
                                          group2_mean=numeric(), 
                                          log10pval=numeric(), 
                                          delta=numeric(), 
                                          condition=character(),
                                          stringsAsFactors=FALSE) 
        }
    }
    DE_table=do.call(rbind,DE_table_list)
    DE_table$condition=substring(DE_table$condition, 2)
    colnames(DE_table)[2:3]=c("target_condition_mean","all_else_mean")
    DE_table$condition=gsub("\\|"," | ",DE_table$condition)
    return(DE_table)
}

#' Create super-cells from Seurat object or a matrix
#' 
#' input a seurat object or single cell matrix, get randomly selected n cells as super-cell seed. For each seed,
#' merge k_merge nearest neighbor cells to the seed cell and calculate average expression.
#'
#' @param sc_object sc matrix with row in cells, col in genes. Or Seurat object, if so Seurat@scale.data will be used for calculation.
#' @param k_filter remove cells that have no mutual nearest neighbor with k=k_filter, default not removing any cells.
#' @param k_merge merge k=k_merge nearest neighbor to create each super-cell.
#' @param n number of super-cell to pick, will be ignored if tag_cell!=NULL.
#' @param tag_cell a vector of cell position for super-cell centers.
#' @param verbose print time consumption to screen or not (0 not print, 1 print important ones, 2 print all).
#' @param sampling_ref a vector group annotation. The random selection of n super-cell centers will try to keep the ratio between group the same.
#' @param seed a integer for set.seed
#' @return an super-cell matrix, each row is a cell, eacho col is a gene. row names is the position of the super-cell center from input (row number in sc matrix, col number in Seurat@scale.data).
#' @export
Super_cell_creation=function(sc_object,k_filter=NULL,k_merge=100,n=5000,tag_cell=NULL,verbose=1,sampling_ref=NULL,seed=123){
    
    set.seed(seed)
    if(typeof(sc_object)=="S4"){
        if(summary(sc_object)[2]=="seurat"){
            pca_out=sc_object@dr$pca@cell.embeddings
            sc_object=t(as.matrix(sc_object@scale.data))
            message("Seurat object detected as input")
        }else{
            stop("invalid sc_object")
        }
    }else{
        message("Assuming input is cell-gene matrix")
        tagMat_filtered = as.matrix( sc_object[,colSums(as.matrix(sc_object) > 0) > ceiling(nrow(sc_object) / 100) ] )
        HVG_idx = order(matrixStats::colSds(tagMat_filtered) / colMeans(tagMat_filtered),decreasing = T)[1:3000]
        pca_out=rsvd::rpca(tagMat_filtered[,HVG_idx],30)$x
    }
    
    if(is.null(sampling_ref)){sampling_ref=nrow(pca_out)}
    
    if(is.null(k_filter)){
        skipfilter=T;k_filter=0
    }else{
        skipfilter=F
    }
    
    ptm <- proc.time()
    
    kmax=max(k_filter,k_merge)
    kmax_plus1=kmax+1
    kfilter_plus1=k_filter+1
    kmerge_plus1=k_merge+1
    
    #mnn construction
    nn_idx=RANN::nn2(pca_out,k=kmax_plus1)$nn.idx
    for(i in 1:nrow(nn_idx)){ #<10s for 20K cells and k=15
        posi=nn_idx[i,] %in% i      
        if(sum(posi)>0){
            posi=which(posi)
            nn_idx[i,posi:kmax]=nn_idx[i,(posi+1):kmax_plus1]
        }
    }
    nn_idx=nn_idx[,-kmax_plus1]
    nn_idx_plusSelf=cbind(1:nrow(nn_idx),nn_idx[,1:k_merge])

    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("MNN done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    
    if(!skipfilter){
        require(igraph)
        edgelist=data.frame(q=rep(1:nrow(pca_out),k_filter),d=(1:nrow(pca_out))[nn_idx[,1:k_filter]])

        ptm <- proc.time()
        colnames(edgelist)=c("from","to")
        edgelist$from=paste0("V",edgelist$from)
        edgelist$to=paste0("V",edgelist$to)

        vertex_mat=data.frame(name=paste0("V",1:nrow(pca_out)),
                              id=rownames(sc_object),
                              stringsAsFactors=F)
        net <- graph_from_data_frame(edgelist, vertices=vertex_mat, directed=T)
        net <- delete_edges(net, which(!which_mutual(net)))
        net <- as.undirected(net,mode="collapse")
        
        #valid targets
        tags=which(degree(net)>0)
        
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        if(verbose){message(paste("network constructed: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    }else{
        tags=1:nrow(pca_out)
    }
    
    if(is.null(tag_cell)){
        #select sample
        if(length(tags)<=n){
            warning(paste0("n bigger than meaningful sample number, n set to ",length(tags)))
            tag_cell=tags
        }else{
            tag_cell=tags[Subsample_by_group(sampling_ref[tags],n)]
        }
    }

    #merging cells to pseudo cells according to target cells
    merging_table=data.frame(center=rep(tag_cell,each=kmerge_plus1),
                             surrounding=as.vector(t(nn_idx_plusSelf[tag_cell,])),
                             stringsAsFactors=F)
    
    #diving merging into multiple chunks
    chunk_cap=30000
    pesudocell_per_chunk=ceiling(chunk_cap/kmerge_plus1)
    row_per_chunk=pesudocell_per_chunk*kmerge_plus1
    number_of_chunk=ceiling(nrow(merging_table)/row_per_chunk)

    ptm <- proc.time()
    chunks=list()
    for(i in 1:number_of_chunk){
        tagrow=((i-1)*row_per_chunk+1):min((i*row_per_chunk),nrow(merging_table))
        chunks[[i]]=as.matrix(sc_object)[merging_table$surrounding[tagrow],]
        chunks[[i]]=cbind(merging_table$center[tagrow],chunks[[i]])
        colnames(chunks[[i]])[1]="cell"
        chunks[[i]]=data.table(chunks[[i]])
        chunks[[i]]=chunks[[i]][,lapply(.SD, mean), by=cell]
        #chunks[[i]]=chunks[[i]][order(chunks[[i]]$cell)]
        rn=chunks[[i]]$cell
        chunks[[i]][,cell:=NULL]
        chunks[[i]]=data.frame(chunks[[i]])
        colnames(chunks[[i]])=colnames(sc_object)
        rownames(chunks[[i]])=rn
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        if(verbose==2){message(paste("chunk ",i,"/",number_of_chunk," finished: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    }
    outmat=do.call(rbind,chunks)
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("merging finished: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    return(outmat)
}
          
#' @export       
plot_table_rank=function(l,tagcol=NULL,topn=10){
    tb=list()
    topElement_list=list()
    for(i in 1:length(l)){
        tb[[i]]=table(l[[i]])
        topElement_list[[i]]=names(tb[[i]])[order(tb[[i]],decreasing = T)[1:topn]]
    }
    topElement=unique(unlist(topElement_list))
    
    if(is.null(tagcol)){
        if(length(l)<=9){
            tagcol=RColorBrewer::brewer.pal(max(3,length(l)),"Set1")
        }else{
            tagcol=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))
            tagcol=tagcol(length(l))
        }
    }
    
    toplot_mat=matrix(0,length(topElement),length(l))
    for(i in 1:length(l)){
        toplot_mat[,i]=tb[[i]][match(topElement,names(tb[[i]]))]/sum(tb[[i]])
    }
    toplot_mat[is.na(toplot_mat)]=0
    
    plot(0,type="n",xlim=c(0,length(topElement)),ylim=c(0,max(toplot_mat)),xlab="",ylab="percent", xaxt='n')
    for(i in 1:length(l)){
        points(toplot_mat[,i],col=tagcol[i],pch=19)
    }
    legend("topright",legend=names(l),bty = "n",lty=0,pch=19,col=tagcol)
    axis(1, at=1:length(topElement), labels=topElement,las=2)
}
