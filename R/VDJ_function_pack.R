Fill_Seurat_DR=function(SeuratOBJ,DRmatrix,DRtype="pca"){
    DR=new("dim.reduction", cell.embeddings = DRmatrix)
    SeuratOBJ@dr[[length(SeuratOBJ@dr)+1]]=DR
    names(SeuratOBJ@dr)[length(SeuratOBJ@dr)]=DRtype
    return(SeuratOBJ)
}

second_to_humanReadableTime=function(t){
  #change second to a vector of hour,min,second
  h=floor(t/3600)
  t=t-h*3600
  m=floor(t/60)
  t=t-m*60
  s=t
  t=c(h,m,s)
  return(t)
}

gp_name_replacing=function(old_group_assignment,old_group_name_to_replace,new_group_name,force_replace=F){
  # replace the group name in old_group_assignment according to old_group_name_to_replace and new_group_name
  # "force_replace = True" allows new_group_name to have the same id as the original group name that is not suppose to be replaced
  all_old_gp_name=unique(old_group_assignment)
  if(!force_replace){
    keeping_group_name=setdiff(all_old_gp_name,old_group_name_to_replace)
    if(any(new_group_name %in% keeping_group_name)){
      stop("there are new group names identical to the original group name that are not supposed to be replaced, set force_replace=T if want to force replace")
    }
  }
  new_group_assignment=old_group_assignment
  tag_unit=which(old_group_assignment %in% old_group_name_to_replace)
  new_group_assignment[tag_unit]=new_group_name[match(old_group_assignment[tag_unit],old_group_name_to_replace)]
  return(new_group_assignment)
}

filter_ident=function(SeuratOBJ,ident_to_rm){
    SeuratOBJ@meta.data$test=rep(0,length(SeuratOBJ@ident))
    SeuratOBJ@meta.data$test[SeuratOBJ@ident %in% ident_to_rm]=1
    if(sum(SeuratOBJ@meta.data$test)>0){
        SeuratOBJ=FilterCells(object = SeuratOBJ, subset.names = "test", high.thresholds = 0.5)
        SeuratOBJ@raw.data=SeuratOBJ@raw.data[,match(colnames(SeuratOBJ@data),colnames(SeuratOBJ@raw.data))]
    }
    return(SeuratOBJ)
}

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

value_spliting=function(in_vector,n=2){
    #split a vector of values by kmean
    tmp_gp=kmeans(in_vector,n)
    gp1=in_vector[tmp_gp$cluster==1]
    gp2=in_vector[tmp_gp$cluster==2]
    threshold=(min(max(gp1),max(gp2))+max(min(gp1),min(gp2)))/2
    return(list(cluster=tmp_gp$cluster,gp_mean=tmp_gp$centers,tmp_gp=threshold))
}

dropcol <- function(df, drop) {
  df <- df [, ! names(df) %in% drop, drop = FALSE]
  return(df)
}

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

write_id_by_group=function(id_list,group_assignment,outprefix){    
    tag_gp=unique(group_assignment)
    tag_gp=tag_gp[gtools::mixedorder(tag_gp)]   
    marker_gene_list=c()
    for(i in 1:length(tag_gp)){
        marker_gene_list[[i]]=id_list[group_assignment==tag_gp[i]]
    }
    #marker_gene_list
    write(paste0(outprefix," seperated by group:"),file=paste0(outprefix,".txt"))
    for(i in 1:length(marker_gene_list)){
        write(paste0(tag_gp[i],":"),file=paste0(outprefix,".txt"), append=T)
        write(marker_gene_list[[i]],file=paste0(outprefix,".txt"), append=T,ncolumns=100000)
        write("",file=paste0(outprefix,".txt"), append=T)
    }
}

split_list_by_group=function(id_list,group_assignment,tag_gp=NULL){  
    if(is.null(tag_gp)){
        tag_gp=unique(group_assignment)
        tag_gp=tag_gp[gtools::mixedorder(tag_gp)] 
    }
    marker_gene_list=c()
    for(i in 1:length(tag_gp)){
        tag=group_assignment==tag_gp[i]
        if(sum(tag)>0){
            marker_gene_list[[i]]=id_list[tag]
            marker_gene_list[[i]]=marker_gene_list[[i]][gtools::mixedorder(marker_gene_list[[i]])]
        }else{
            marker_gene_list[[i]]=c()
        }
    }
    names(marker_gene_list)=tag_gp
    return(marker_gene_list)
}

lsnofun <- function(name = parent.frame()) {
    obj <- ls(name = name)
    obj[!sapply(obj, function(x) is.function(get(x)))]
}

DEG_wilcox=function(group1,group2,delta_threshold=NULL,p_threshold=NULL,p_adjust_method='bonferroni'){
  #DEG for UMI or TPM (value >0)
  if(sum(colnames(group1)!=colnames(group2))>0){
    warning("gene names are different between group1 and group2")
  }
  p=rep(1,ncol(group1))
  g1_mean=colMeans(group1)
  g2_mean=colMeans(group2)
  mean_dif=g1_mean-g2_mean
  tag=which(mean_dif!=0)
  p[tag]=sapply(tag,function(x){wilcox.test(group1[,x], group2[,x],exact=F)$p.value})
  #for fold change 0 sum not allowed

  p=p.adjust(p,method=p_adjust_method)
  #check threshold
  if(!is.null(p_threshold) & !is.null(delta_threshold)){
    tag=p<p_threshold & abs(mean_dif) > abs(delta_threshold)
    o=data.frame(group1_mean=g1_mean,group2_mean=g2_mean,log10pval=log10(p),delta=mean_dif,DE=tag,stringsAsFactors=F)
  }else{
    o=data.frame(group1_mean=g1_mean,group2_mean=g2_mean,log10pval=log10(p),delta=mean_dif,stringsAsFactors=F)
  }
  rownames(o)=colnames(group1)
  return(o)
}                
                
                
DEG_wilcox_UMI=function(group1,group2,p_threshold=NULL,log2fold_threshold=NULL,p_adjust_method='fdr'){ #it was bonferroni
  #DEG for UMI or TPM (value >0)
  if(sum(colnames(group1)!=colnames(group2))>0){
    warning("gene names are different between group1 and group2")
  }
  p=rep(1,ncol(group1))
  foldchange=rep(1,ncol(group1))
  g1_mean=colMeans(group1)
  g2_mean=colMeans(group2)
  g1_pos=g1_mean>0
  g2_pos=g2_mean>0
  tag=which(g1_pos | g2_pos)
  p[tag]=sapply(tag,function(x){wilcox.test(group1[,x], group2[,x],exact=F)$p.value})
  #for fold change 0 sum not allowed
  tag=which(g1_pos | g2_pos)
  foldchange[tag]=sapply(tag,function(x){log2((g1_mean[x]+1)/(g2_mean[x]+1))})
  foldchange[which(g1_pos>g2_pos)]=Inf
  foldchange[which(g1_pos<g2_pos)]=-Inf
  
  p=p.adjust(p,method=p_adjust_method)
  #check threshold
  if(!is.null(p_threshold) & !is.null(log2fold_threshold)){
    tag=p<p_threshold & abs(foldchange) > abs(log2fold_threshold)
    o=data.frame(group1_mean=g1_mean,group2_mean=g2_mean,log10pval=log10(p),log2fold=foldchange,DE=tag,stringsAsFactors=F)
  }else{
    o=data.frame(group1_mean=g1_mean,group2_mean=g2_mean,log10pval=log10(p),log2fold=foldchange,stringsAsFactors=F)
  }
  rownames(o)=colnames(group1)
  return(o)
}
                
matrix_Aggregate=function(mat,rowby,colby,function_to_use){
  #function_to_use: mean, median, max, min, etc
  dt=data.table(cbind(rowby,mat))
  dt=dt[,lapply(.SD, function(x){function_to_use(x,na.rm = T)}),by=rowby]
  dt=dt[,rowby:=NULL]
  dt=t(dt)
  dt=data.table(cbind(colby,dt))
  dt=dt[,lapply(.SD, function(x){function_to_use(x,na.rm = T)}),by=colby]
  dt=dt[,colby:=NULL]
  dt=as.matrix(t(dt))
  return(dt)
}
                
cross_combining=function(string1,string2){
    o=rep("",length(string1)*length(string2))
    k=0
    for(i in 1:length(string1)){
        for(j in 1:length(string2)){
            k=k+1
            o[k]=paste0(string1[i],"_",string2[j])
        } 
    }
    return(o)
}
                
Subsample_by_group=function(assign_vector,tag_total_num){
    ori_tot=length(assign_vector)
    if(tag_total_num>=ori_tot){
        message("no need to down sample")
        return(1:length(assign_vector))
    }else{
        types=unique(assign_vector)
        tags=c()
        for(i in 1:length(types)){
            in_gp=which(assign_vector==types[i])
            tags=c(tags,sample(in_gp,round(length(in_gp)/ori_tot*tag_total_num))) 
        }
        return(tags)
    }
}

Subsample_by_group_and_importance=function(assign_vector,importance_score,tag_total_num){
    #same as subsample by group, but select the top n important target for each group
    ori_tot=length(assign_vector)
    if(tag_total_num>=ori_tot){
        message("no need to down sample")
        return(1:length(assign_vector))
    }else{
        types=unique(assign_vector)
        tags=c()
        for(i in 1:length(types)){
            in_gp=which(assign_vector==types[i])
            tags=c(tags,in_gp[order(importance_score[in_gp],decreasing=T)[1:(round(length(in_gp)/ori_tot*tag_total_num))]]) 
        }
        return(tags)
    }
}                


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
                
VJC_stats=function(position_index_matrix,V_id,J_id,C_id,batch_info){
    tag=which(position_index_matrix$v_gene==V_id & position_index_matrix$j_gene==J_id & position_index_matrix$c_gene==C_id)
    #tag=tag[!duplicated(position_index_matrix$barcode_idx[tag])]
    #message(paste0("total number of cells of given VJC:",length(unique(position_index_matrix$barcode_idx[tag]))))

    a1=table(position_index_matrix$celltype[tag],useNA="i")
    a2=table(position_index_matrix$celltype[!duplicated(position_index_matrix$barcode)],useNA="i")
    a1=a1[match(names(a2),names(a1))]
    a1[is.na(a1)]=0
    a1=cbind(count=a1,percentage=paste0(round(a1/a2*100,2),"%"))
    rownames(a1)=names(a2)
    rownames(a1)[is.na(rownames(a1))]="NA"
    write.xlsx(data.frame(a1), file=paste0(V_id,"_",J_id,"_",C_id,"_stats.xlsx"), sheetName="celltype", append=F, row.names=T,showNA = T)

    a1=table(position_index_matrix$batch_id[tag],useNA="i")
    a2=table(position_index_matrix$batch_id[!duplicated(position_index_matrix$barcode)],useNA="i")
    a1=a1[match(names(a2),names(a1))]
    a1[is.na(a1)]=0
    a1=cbind(count=a1,percentage=a1/a2) #paste0(round(a1/a2*100,2),"%")
    rownames(a1)=names(a2)

    batch_info=batch_info[match(rownames(a1),batch_info$batch),]
    a1=data.table(cbind(a1,batch_info))
    b1=a1[,.(count=sum(count),percentage=mean(percentage)),by=patient_id]
    c1=a1[,.(count=sum(count),percentage=mean(percentage)),by=tissue_type]
    d1=a1[,.(count=sum(count),percentage=mean(percentage)),by=disease_status]

    a1_w=data.frame(count=a1$count,percentage=paste0(round(a1$percentage*100,2),"%"))
    rownames(a1_w)=a1$batch
    b1_w=data.frame(count=b1$count,percentage=paste0(round(b1$percentage*100,2),"%"))
    rownames(b1_w)=b1$patient_id
    c1_w=data.frame(count=c1$count,percentage=paste0(round(c1$percentage*100,2),"%"))
    rownames(c1_w)=c1$tissue_type
    d1_w=data.frame(count=d1$count,percentage=paste0(round(d1$percentage*100,2),"%"))
    rownames(d1_w)=d1$disease_status

    write.xlsx(data.frame(a1_w), file=paste0(V_id,"_",J_id,"_",C_id,"_stats.xlsx"), sheetName="batch_id", append=T, row.names=T,showNA = T)
    write.xlsx(data.frame(b1_w), file=paste0(V_id,"_",J_id,"_",C_id,"_stats.xlsx"), sheetName="patient_id", append=T, row.names=T,showNA = T)
    write.xlsx(data.frame(c1_w), file=paste0(V_id,"_",J_id,"_",C_id,"_stats.xlsx"), sheetName="tissue_type", append=T, row.names=T,showNA = T)
    write.xlsx(data.frame(d1_w), file=paste0(V_id,"_",J_id,"_",C_id,"_stats.xlsx"), sheetName="disease_status", append=T, row.names=T,showNA = T)
}
                
plot_usage=function(position_index_matrix,tag_col,out_name,out_height,out_width,lhei,lwid,margins,notecex,Colv=T,Rowv=T,colorby=c("row","col","raw")){
    col2=rev(RColorBrewer::brewer.pal(9,"RdBu"))
        
    pdf(out_name,height=out_height,width=out_width)
    tmp=table(position_index_matrix$celltype,position_index_matrix[,tag_col])
    if(colorby=="row"){
        toplot=t(apply(tmp,1,function(x){x/sum(x)}))
    }else if(colorby=="col"){
        toplot=apply(tmp,2,function(x){x/sum(x)})
    }else if(colorby=="raw"){
        toplot=tmp
    }else{
        stop("colorby must be row / col / raw")
    }
    heatmap.2(toplot,col=colorRampPalette(col2),
              trace="none",density.info="none",dendrogram='none',Colv=Colv,Rowv=Rowv,notecol="black",key=T,lhei=lhei,lwid=lwid,
              main=paste0("Cell type by ",tag_col),key.title="",key.xlab="",margins=margins,cellnote=tmp,notecex=notecex)
    
    tmp=table(position_index_matrix$batch_id,position_index_matrix[,tag_col])
    if(colorby=="row"){
        toplot=t(apply(tmp,1,function(x){x/sum(x)}))
    }else if(colorby=="col"){
        toplot=apply(tmp,2,function(x){x/sum(x)})
    }else if(colorby=="raw"){
        toplot=tmp
    }else{
        stop("colorby must be row / col / raw")
    }
    heatmap.2(toplot,col=colorRampPalette(col2),
              trace="none",density.info="none",dendrogram='none',Colv=Colv,Rowv=Rowv,notecol="black",key=T,lhei=lhei,lwid=lwid,
              main=paste0("Sample ID by ",tag_col),key.title="",key.xlab="",margins=margins,cellnote=tmp,notecex=notecex)

    tmp=table(position_index_matrix$tissue_type,position_index_matrix[,tag_col])
    if(colorby=="row"){
        toplot=t(apply(tmp,1,function(x){x/sum(x)}))
    }else if(colorby=="col"){
        toplot=apply(tmp,2,function(x){x/sum(x)})
    }else if(colorby=="raw"){
        toplot=tmp
    }else{
        stop("colorby must be row / col / raw")
    }
    heatmap.2(toplot,col=colorRampPalette(col2),
              trace="none",density.info="none",dendrogram='none',Colv=Colv,Rowv=Rowv,notecol="black",key=T,lhei=lhei,lwid=lwid,
              main=paste0("Tissue type by ",tag_col),key.title="",key.xlab="",margins=margins,cellnote=tmp,notecex=notecex)

    tmp=table(position_index_matrix$patient_id,position_index_matrix[,tag_col])
    if(colorby=="row"){
        toplot=t(apply(tmp,1,function(x){x/sum(x)}))
    }else if(colorby=="col"){
        toplot=apply(tmp,2,function(x){x/sum(x)})
    }else if(colorby=="raw"){
        toplot=tmp
    }else{
        stop("colorby must be row / col / raw")
    }
    heatmap.2(toplot,col=colorRampPalette(col2),
              trace="none",density.info="none",dendrogram='none',Colv=Colv,Rowv=Rowv,notecol="black",key=T,lhei=lhei,lwid=lwid,
              main=paste0("Patient ID by ",tag_col),key.title="",key.xlab="",margins=margins,cellnote=tmp,notecex=notecex)

    tmp=table(position_index_matrix$disease,position_index_matrix[,tag_col])
    if(colorby=="row"){
        toplot=t(apply(tmp,1,function(x){x/sum(x)}))
    }else if(colorby=="col"){
        toplot=apply(tmp,2,function(x){x/sum(x)})
    }else if(colorby=="raw"){
        toplot=tmp
    }else{
        stop("colorby must be row / col / raw")
    }
    heatmap.2(toplot,col=colorRampPalette(col2),
              trace="none",density.info="none",dendrogram='none',Colv=Colv,Rowv=Rowv,notecol="black",key=T,lhei=lhei,lwid=lwid,
              main=paste0("Disease status by ",tag_col),key.title="",key.xlab="",margins=margins,cellnote=tmp,notecex=notecex)
    dev.off()
}

spreading_cal=function(tag_patient,tag_tissue,position_index_matrix){
    #clonetype that expand from one tissue to the other
    expansion_info=list()
    no_expansion_clones=c()
    locally_expanded_clones=c()
    cross_tissue_clones=c()
    for(i in 1:length(tag_patient)){
        tag=position_index_matrix$patient_id==tag_patient[i] & position_index_matrix$tissue_type %in% tag_tissue
        batch_clontypeid_tbl=table(position_index_matrix$batch[tag],position_index_matrix$clone_id[tag])
        total_cell_perbatch=rowSums(batch_clontypeid_tbl)
        batch_clontypeid_tbl_noexp=batch_clontypeid_tbl[,colSums(batch_clontypeid_tbl)==1,drop=FALSE]
        no_expansion_clones=c(no_expansion_clones,colnames(batch_clontypeid_tbl)[colSums(batch_clontypeid_tbl)==1])
        batch_clontypeid_tbl=batch_clontypeid_tbl[,colSums(batch_clontypeid_tbl)>1,drop=FALSE]    #update
        batch_clontypeid_tbl_explocally=batch_clontypeid_tbl[,colSums(batch_clontypeid_tbl>0)==1,drop=FALSE]
        locally_expanded_clones=c(locally_expanded_clones,colnames(batch_clontypeid_tbl)[colSums(batch_clontypeid_tbl>0)==1])
        batch_clontypeid_tbl=batch_clontypeid_tbl[,colSums(batch_clontypeid_tbl>0)>1,drop=FALSE]    #update
        cross_tissue_clones=c(cross_tissue_clones,colnames(batch_clontypeid_tbl))
        expansion_info[[i]]=data.frame(batch_id=rownames(batch_clontypeid_tbl_noexp),
                                       no_expansion=rowSums(batch_clontypeid_tbl_noexp),
                                       expanded_locally=rowSums(batch_clontypeid_tbl_explocally),
                                       expanded_cross_tissue=rowSums(batch_clontypeid_tbl))
    }
    expansion_info_tbl=data.table(do.call(rbind,expansion_info))
    expansion_info_tbl_long= melt(expansion_info_tbl, id.vars = c("batch_id"),
                                  measure.vars = c("no_expansion", "expanded_locally", "expanded_cross_tissue"))
    expansion_info_tbl_long$patient=gsub("_.*","",unique(expansion_info_tbl_long$batch_id))
    expansion_info_tbl_long$tissue=gsub("^.*_","",unique(expansion_info_tbl_long$batch_id))
    
    return(list(tbl_wide=expansion_info_tbl,
                tbl_long=expansion_info_tbl_long,
                no_expansion_clones=no_expansion_clones,
                locally_expanded_clones=locally_expanded_clones,
                cross_tissue_clones=cross_tissue_clones))
}
        
expansion_per_condition=function(target_condition,position_index_matrix,plot_width,plot_height,TorB="Tcell"){
    tmp=table(eval(parse(text = paste0("position_index_matrix$",target_condition))),position_index_matrix$clone_id)
    to_plot=data.table(target_condition=rownames(tmp),
                       no_expansion=0,
                       low_expansion=0,
                       high_expansion=0)
    cellnum_per_row=rowSums(tmp)
    
    for(i in 1:nrow(tmp)){
        to_plot$no_expansion[i]=sum(tmp[i,]==1)
        to_plot$low_expansion[i]=sum(tmp[i,tmp[i,]>1 & (tmp[i,]<=cellnum_per_row[i]/100 | tmp[i,]<=2)])
        to_plot$high_expansion[i]=sum(tmp[i,tmp[i,]>2 & tmp[i,]>cellnum_per_row[i]/100])
    }
    to_plot_long= melt(to_plot, id.vars = "target_condition",
                       measure.vars = c("no_expansion", "low_expansion", "high_expansion"),)
    
    g=ggplot() +
        geom_bar(aes(y = value, x = target_condition, fill = variable), data = to_plot_long,stat="identity",position = "fill") +
        scale_y_continuous(labels = percent_format()) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        labs(fill='expansion level',x=target_condition,y="percentage of total cell number")
    pdf(paste0(TorB,"_expansion_lvl_by_",target_condition,".pdf"),width=plot_width,height=plot_height)
        print(g)
    dev.off()
    
    return(list(tbl_wide=to_plot,tbl_long=to_plot_long))
}
        
cell_type_distribution_in_spreading=function(position_index_matrix,expansion_type_list,tag_patient,tag_tissue,celltype_assignment_all){
    tag_batch=cross_combining(tag_patient,tag_tissue)
    
    tmp=list()
    for(i in 1:length(tag_batch)){    
        tag=position_index_matrix$batch_id==tag_batch[i]
        expansion_type_cell_type_table=table(expansion_type_list[tag],position_index_matrix$celltype[tag],useNA = "ifany") #
        colnames(expansion_type_cell_type_table)[is.na(colnames(expansion_type_cell_type_table))]="unknown"
        if(sum(expansion_type_cell_type_table[,-ncol(expansion_type_cell_type_table)])==0){
            next
        }
        #expansion_type_cell_type_table[,-ncol(expansion_type_cell_type_table)]=expansion_type_cell_type_table[,-ncol(expansion_type_cell_type_table)]/sum(expansion_type_cell_type_table[,-ncol(expansion_type_cell_type_table)])
        #expansion_type_cell_type_table[,"unknown"]=expansion_type_cell_type_table[,"unknown"]/sum(expansion_type_cell_type_table)
        expansion_type_cell_type_table=expansion_type_cell_type_table[,-ncol(expansion_type_cell_type_table)]
        expansion_type_cell_type_table=expansion_type_cell_type_table/sum(expansion_type_cell_type_table)
        expansion_type_cell_type_df=as.data.frame.matrix(expansion_type_cell_type_table)
        expansion_type_cell_type_df$expansion_type=rownames(expansion_type_cell_type_df)
        expansion_type_cell_type_df=data.table(expansion_type_cell_type_df)
        expansion_type_cell_type_df= melt(expansion_type_cell_type_df, id.vars = c("expansion_type"),
                                          measure.vars = setdiff(colnames(expansion_type_cell_type_df),"expansion_type"),
                                          variable.name = "cell_type", value.name = "percent_of_patient")
        expansion_type_cell_type_df$batch=tag_batch[i]
        tmp[[i]]=expansion_type_cell_type_df
    }
    toplot=do.call(rbind,tmp)
    toplot$patient=gsub("_.*","",toplot$batch)
    toplot$tissue=gsub("^.*_","",toplot$batch)
    
    Group_type=unique(celltype_assignment_all)
    Group_type=Group_type[gtools::mixedorder(Group_type)]    
    present_celltype=unique(toplot$cell_type)    
    coltouse=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(length(Group_type))
    coltouse=coltouse[Group_type %in% present_celltype]
    
    g=ggplot() +
            geom_bar(aes(y = percent_of_patient, x = patient, fill = factor(cell_type,levels=Group_type)), data = toplot,stat="identity",position = "stack") +
            scale_y_continuous(labels = percent_format()) +
            theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
            facet_wrap( ~ tissue + expansion_type,scales="free") + 
            ggtitle(paste(tag_tissue,collapse=" & ")) + labs(y = "pecentage of cells in given library",fill='cell type') +
            scale_fill_manual("legend", values = coltouse)
    return(list(toplot=toplot,g=g))               
}
        
umap_seurat=function(seurat_in,pca_dim=25){
    umap_out=umap::umap(seurat_in@dr$pca@cell.embeddings[,1:pca_dim])
    seurat_in@dr$umap=seurat_in@dr$tsne
    seurat_in@dr$umap@key='umap_'
    seurat_in@dr$umap@cell.embeddings=umap_out$layout
    return(seurat_in)
}
        
readDICE=function(DICE_path="/home/ahe/Analysis/201801_JohnVDJ/data/DICE_database"){
    #DICE_path="/home/ahe/Analysis/201801_JohnVDJ/data/DICE_database"
    file_names=list.files(path=DICE_path,pattern="_TPM.csv",full.names=F)
    celltypes=gsub("_TPM.csv","",file_names)
    DICE_tpm=list()
    for(i in 1:length(file_names)){
        DICE_tpm[[i]]=fread(paste0(DICE_path,"/",file_names[i]),sep=",",data.table=F)
    }
    rname=gsub(";.*","",DICE_tpm[[1]][,3])
    k=1
    dup=duplicated(rname)
    while(any(dup)){
        rname[dup]=paste0(rname[dup],"_",k)
        k=k+1
        dup=duplicated(rname)
    }
    ave_exp=matrix(0,nrow(DICE_tpm[[i]]),length(celltypes))
    rownames(ave_exp)=rname
    colnames(ave_exp)=celltypes
    for(i in 1:length(file_names)){
        DICE_tpm[[i]]=DICE_tpm[[i]][,-c(1,2,3,4)]
        rownames(DICE_tpm[[i]])=rname
        ave_exp[,i]=rowMeans(DICE_tpm[[i]][,-c(1,2,3,4)])
    }
    names(DICE_tpm)=celltypes
    return(list(ave_exp=ave_exp,exp_list=DICE_tpm))
}

cb=function(in1,in2){
    out=list()
    k=1
    for(i in 1:length(in1)){
        for(j in 1:length(in2)){
            out[[k]]=c(unlist(in1[i]),unlist(in2[j]))
            k=k+1
        }
    }
    return(out)
}

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
        
Sampling_merging_seurat=function(tagSeurat,k_filter=15,k_merge=100,n=5000,tag_cell=NULL,verbose=1,sampling_ref=rep(1,ncol(tagSeurat@data))){
    #input a seurat object (tagSeurat), get randomly selected n cells as super cell seed. For each seed, merge k_merge nearest neighbor cells to the seed cell and calculate average expression.
    # tagSeurat@scale.data were used for calculation.
    #setup

    if(is.null(k_filter)){
        skipfilter=T;k_filter=0
    }else{
        skipfilter=F
    }
    
    ptm <- proc.time()
    pca_out=tagSeurat@dr$pca@cell.embeddings
    kmax=max(k_filter,k_merge)
    kmax_plus1=kmax+1
    kfilter_plus1=k_filter+1
    kmerge_plus1=k_merge+1
    
    #mnn construction
    nn_idx=nn2(pca_out,k=kmax_plus1)$nn.idx
    for(i in 1:nrow(nn_idx)){ #<10s for 20K cells and k=15
        posi=nn_idx[i,] %in% i      
        if(sum(posi)>0){
            posi=which(posi)
            nn_idx[i,posi:kmax]=nn_idx[i,(posi+1):kmax_plus1]
        }
    }
    nn_idx=nn_idx[,-kmax_plus1]
    nn_idx_plusSelf=cbind(1:nrow(nn_idx),nn_idx[,1:k_merge])

    if(!skipfilter){
        edgelist=data.frame(q=rep(1:nrow(pca_out),k_filter),d=(1:nrow(pca_out))[nn_idx[,1:k_filter]])
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        if(verbose){message(paste("MNN done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}

        ptm <- proc.time()
        colnames(edgelist)=c("from","to")
        edgelist$from=paste0("V",edgelist$from)
        edgelist$to=paste0("V",edgelist$to)

        vertex_mat=data.frame(name=paste0("V",1:nrow(pca_out)),
                              id=rownames(tagSeurat@meta.data),
                              #sample_gp=tagSeurat@meta.data$sample_gp,
                              #patient_assignment=tagSeurat@meta.data$patient_assignment,
                              #tissue_assignment=tagSeurat@meta.data$tissue_assignment,
                              #disease_assignment=tagSeurat@meta.data$disease_assignment,
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
        chunks[[i]]=t(as.matrix(tagSeurat@scale.data))[merging_table$surrounding[tagrow],]
        chunks[[i]]=cbind(merging_table$center[tagrow],chunks[[i]])
        colnames(chunks[[i]])[1]="cell"
        chunks[[i]]=data.table(chunks[[i]])
        chunks[[i]]=chunks[[i]][,lapply(.SD, mean), by=cell]
        #chunks[[i]]=chunks[[i]][order(chunks[[i]]$cell)]
        rn=chunks[[i]]$cell
        chunks[[i]][,cell:=NULL]
        chunks[[i]]=data.frame(chunks[[i]])
        colnames(chunks[[i]])=rownames(as.matrix(tagSeurat@scale.data))
        rownames(chunks[[i]])=rn
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        if(verbose==2){message(paste("chunk ",i,"/",number_of_chunk," finished: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    }
    tagMat=do.call(rbind,chunks)
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("merging finished: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    return(tagMat)
}
        
Sampling_merging=function(tagMat,k_filter=15,k_merge=100,n=5000,tag_cell=NULL,verbose=1){
    #input a cell-gene matrix (tagMat), get randomly selected n cells as super cell seed. For each seed, merge k_merge nearest neighbor cells to the seed cell and calculate average expression. 
    tagMat_filtered=as.matrix(tagMat[,colSums(as.matrix(tagMat)>0)>ceiling(nrow(tagMat)/100)])
    
    #setup
    ptm <- proc.time()
    #variable gene
    HVG_idx=order(matrixStats::colSds(tagMat_filtered)/colMeans(tagMat_filtered),decreasing = T)[1:3000]
    
    pca_out=rsvd::rpca(tagMat_filtered[,HVG_idx],30)$x
    kmax=max(k_filter,k_merge)
    kmax_plus1=kmax+1
    kfilter_plus1=k_filter+1
    kmerge_plus1=k_merge+1
    
    #mnn construction
    nn_idx=nn2(pca_out,k=kmax_plus1)$nn.idx
    for(i in 1:nrow(nn_idx)){ #<10s for 20K cells and k=15
        posi=nn_idx[i,] %in% i      
        if(sum(posi)>0){
            posi=which(posi)
            nn_idx[i,posi:kmax]=nn_idx[i,(posi+1):kmax_plus1]
        }
    }
    nn_idx=nn_idx[,-kmax_plus1]
    nn_idx_plusSelf=cbind(1:nrow(nn_idx),nn_idx[,1:k_merge])

    edgelist=data.frame(q=rep(1:nrow(pca_out),k_filter),d=(1:nrow(pca_out))[nn_idx[,1:k_filter]])
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("MNN done: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}

    ptm <- proc.time()
    colnames(edgelist)=c("from","to")
    edgelist$from=paste0("V",edgelist$from)
    edgelist$to=paste0("V",edgelist$to)

    vertex_mat=data.frame(name=paste0("V",1:nrow(pca_out)),
                          id=rownames(tagMat),
                          stringsAsFactors=F)
    net <- graph_from_data_frame(edgelist, vertices=vertex_mat, directed=T)
    net <- delete_edges(net, which(!which_mutual(net)))
    net <- as.undirected(net,mode="collapse")
    #net <- simplify(net)
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("network constructed: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    
    if(is.null(tag_cell)){
        #valid targets
        tags=which(degree(net)>0)

        #select sample
        if(length(tags)<=n){
            warning(paste0("n bigger than meaningful sample number, n set to ",length(tags)))
            tag_cell=tags
        }else{
            tag_cell=sample(tags,n)
        }
    }

    #merging cells to pseudo cells according to target cells
    merging_table=data.frame(center=rep(tag_cell,each=kmerge_plus1),
                             surrounding=as.vector(t(nn_idx_plusSelf[tag_cell,])),
                             stringsAsFactors=F)
    
    #diving merging into multiple chunks
    chunk_cap=40000
    pesudocell_per_chunk=ceiling(chunk_cap/kmerge_plus1)
    row_per_chunk=pesudocell_per_chunk*kmerge_plus1
    number_of_chunk=ceiling(nrow(merging_table)/row_per_chunk)

    ptm <- proc.time()
    chunks=list()
    for(i in 1:number_of_chunk){
        tagrow=((i-1)*row_per_chunk+1):min((i*row_per_chunk),nrow(merging_table))
        chunks[[i]]=as.matrix(tagMat)[merging_table$surrounding[tagrow],]
        chunks[[i]]=cbind(merging_table$center[tagrow],chunks[[i]])
        colnames(chunks[[i]])[1]="cell"
        chunks[[i]]=data.table(chunks[[i]])
        chunks[[i]]=chunks[[i]][,lapply(.SD, mean), by=cell]
        #chunks[[i]]=chunks[[i]][order(chunks[[i]]$cell)]
        rn=chunks[[i]]$cell
        chunks[[i]][,cell:=NULL]
        chunks[[i]]=data.frame(chunks[[i]])
        colnames(chunks[[i]])=colnames(as.matrix(tagMat))
        rownames(chunks[[i]])=rn
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        if(verbose==2){message(paste("chunk ",i,"/",number_of_chunk," finished: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    }
    tagMat=do.call(rbind,chunks)
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    if(verbose){message(paste("merging finished: time consumed:",t[1],"hr",t[2],"min",t[3],"s"))}
    return(tagMat)
}

bulk_sample_comparison=function(bulk_df,sc_mat,cor_threshold=NULL,scale_cor_2nd=T,base_coefficient=NULL){
    #both of the expression_matrix, row as sample, col as feature
    #colnames identical
    mean_bulk=colMeans(bulk_df)
    mean_extracted_bulk=scale(bulk_df,scale = F,center = T)

    cor_raw=cor(t(sc_mat),mean_bulk)
    
    if(is.null(cor_threshold)){
        options(repr.plot.width=6, repr.plot.height=3)
        plot(density(cor_raw))
        abline(v=seq(-1,1,0.1),lty=2)
        return(NULL)
    }else{
        mean_SC_for_Extraction=colMeans(sc_mat[cor_raw[,1]>cor_threshold,])
        mean_extracted_SC=t(apply(sc_mat,1,function(x){x-mean_SC_for_Extraction}))

        cor_2nd=cor(t(mean_extracted_SC),t(mean_extracted_bulk))
        
        if(scale_cor_2nd){
            if(is.null(base_coefficient)){
                base_coefficient=cor_raw-cor_threshold
            }
            base_coefficient[base_coefficient<0]=0

            for(i in 1:ncol(cor_2nd)){
                cor_2nd[,i]=cor_2nd[,i]*base_coefficient
            }
        }

        return(list(cor_raw=cor_raw,
                    cor_2nd=cor_2nd,
                    base_coefficient=base_coefficient))
    }
}
        
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
    return(toplot_mat)
}