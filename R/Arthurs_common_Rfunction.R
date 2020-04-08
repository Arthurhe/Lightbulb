#' @export
second_to_humanReadableTime=function(t){
  #change second to a vector of hour,min,second
  h=floor(t/3600)
  t=t-h*3600
  m=floor(t/60)
  t=t-m*60
  s=round(t,2)
  t=c(h,m,s)
  return(t)
}

#' @export
bedtools_intersect=function(bed1,bed2,overlap=T,count=F,maxgap=0,minoverlap=0,ignore_strand=T,strand_col=6){
    suppressMessages(require(GenomicRanges))
    if(ignore_strand){
        bed1 <- GRanges(seqnames = bed1[,1],
                        ranges = IRanges(start = bed1[,2],
                                         end = bed1[,3]))
        bed2 <- GRanges(seqnames = bed2[,1],
                        ranges = IRanges(start = bed2[,2],
                                         end = bed2[,3]))        
    }else{
        bed1 <- GRanges(seqnames = bed1[,1],
                        ranges = IRanges(start = bed1[,2],
                                         end = bed1[,3]),
                        strand=bed1[,strand_col])
        bed2 <- GRanges(seqnames = bed2[,1],
                        ranges = IRanges(start = bed2[,2],
                                         end = bed2[,3]),
                        strand=bed2[,strand_col])          
    }
    if(overlap){
        overlap=findOverlaps(bed1,bed2,ignore.strand = ignore_strand,maxgap=maxgap,minoverlap=minoverlap)
        overlap=data.frame(overlap)
        colnames(overlap)=c("bed1_idx","bed2_idx")
    }
    if(count){
        count=countOverlaps(bed1,bed2,ignore.strand = ignore_strand,maxgap=maxgap,minoverlap=minoverlap)
    }
    return(list(overlap=overlap,count=count))
}

#' @export
gp_name_replacing=function(old_group_assignment,old_group_name_to_replace,new_group_name,force_replace=F){
  # replace the group name in old_group_assignment according to old_group_name_to_replace and new_group_name
  # "force_replace = True" allows new_group_name to have the same id as the original group name that is not suppose to be replaced
    if(length(old_group_name_to_replace)!=length(new_group_name)){stop("gp names length not identical")}
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

#' @export
split_list_by_group=function(id_list,group_assignment,tag_gp=NULL){  
    if(is.null(tag_gp)){
        tag_gp=unique(group_assignment)
        tag_gp=tag_gp[gtools::mixedorder(tag_gp)] 
    }
    marker_gene_list=list()
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

#' @export
lsnofun <- function(name = parent.frame()) {
    obj <- ls(name = name)
    obj[!sapply(obj, function(x) is.function(get(x)))]
}
                
#' @export                 
Subsample_by_group_BreakEven=function(assign_vector,tag_total_num){
    ori_tot=length(assign_vector)
    if(tag_total_num>=ori_tot){
        message("no need to down sample")
        return(1:length(assign_vector))
    }else{
        types=unique(assign_vector)
        numpertype=ceiling(tag_total_num/length(types))
        tags=c()
        for(i in 1:length(types)){
            in_gp=which(assign_vector==types[i])
            if(length(in_gp)<=numpertype){
                tags=c(tags,in_gp) 
            }else{
                tags=c(tags,sample(in_gp,numpertype)) 
            }
        }
        return(tags)
    }
}
                
#' @export                 
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
                
#' @export 
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
                
#' @export                 
which.colmax=function(input_mat){
    return(apply(input_mat,2,function(x){which.max(x)}))
}
                
#' @export 
matrix_Aggregate=function(mat,rowby,colby,function_to_use){
  #function_to_use: mean, median, max, min, etc
  rowfactor=factor(rowby)  
  colfactor=factor(colby)
  rowby=as.numeric(rowfactor)
  colby=as.numeric(colfactor)
  dt=data.table(cbind(rowby,mat))
  dt=dt[,lapply(.SD, function(x){function_to_use(x,na.rm = T)}),by=rowby]
  dt=dt[order(dt$rowby),]
  rname=dt$rowby  
  dt=dt[,rowby:=NULL]
  dt=t(dt)
  dt=data.table(cbind(colby,dt))
  dt=dt[,lapply(.SD, function(x){function_to_use(x,na.rm = T)}),by=colby]
  dt=dt[order(dt$colby),]
  cname=dt$colby
  dt=dt[,colby:=NULL]
  dt=as.matrix(t(dt))
  rownames(dt)=levels(rowfactor)[rname]
  colnames(dt)=levels(colfactor)[cname]  
  return(dt)
}
                
#' @export                 
firstup <- function(x) {
    x=tolower(x)
   substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    return(x)
}
                
#' @export 
colSds=function(x){
    return(matrixStats::colSds(as.matrix(x)))
}
                
#' @export                 
rowSds=function(x){
    return(matrixStats::rowSds(as.matrix(x)))
}
                
#' @export               
colMaxs=function(x){
    return(matrixStats::colMaxs(as.matrix(x)))
}
                
#' @export                 
rowMaxs=function(x){
    return(matrixStats::rowMaxs(as.matrix(x)))
}
                
#' @export                 
colMins=function(x){
    return(matrixStats::colMins(as.matrix(x)))
}
                
#' @export                 
rowMins=function(x){
    return(matrixStats::rowMins(as.matrix(x)))
}
                
#' @export                 
reassign_cluster=function(SNN,cluster_annotation,unwannted_cluster){ #SNN can be any similarity matrix that the higher the score the higher similarity (lower distance)
    unassigned_cells=which(cluster_annotation %in% unwannted_cluster)
    assigned_cells=which(!cluster_annotation %in% unwannted_cluster)
    loop_counter=1
    while(length(unassigned_cells)>0 & loop_counter<100){
        sub_SNN=SNN[unassigned_cells,assigned_cells]
        new_assignment=cluster_annotation[unassigned_cells]
        
        if(length(unassigned_cells)==1){
            if(max(sub_SNN)>0){
                cluster_annotation[unassigned_cells]=cluster_annotation[assigned_cells][which.max(sub_SNN)]
                break()
            }else{
                cluster_annotation[unassigned_cells]=NA
            }
        }
        
        new_assigned_cells=which(rowMaxs(sub_SNN)>0)
        if(length(new_assigned_cells)==0){
            cluster_annotation[unassigned_cells]=NA
            break()
        }else{
            for(i in 1:length(new_assigned_cells)){
                new_assignment[new_assigned_cells[i]]=cluster_annotation[assigned_cells][which.max(sub_SNN[new_assigned_cells[i],])]
            }
        }

        #put new assignment to annotation
        cluster_annotation[unassigned_cells]=new_assignment
        
        #renew unassigned cell list
        assigned_cells=c(assigned_cells,unassigned_cells[new_assigned_cells])
        unassigned_cells=unassigned_cells[-new_assigned_cells]
        loop_counter=loop_counter+1
    }
    return(cluster_annotation)
}

           
#' @export
DEG=function(group1,group2,delta_threshold=NULL,p_threshold=NULL,p_adjust_method='bonferroni',test.method=c("wilcox","t")){
  #DEG for UMI or TPM (value >0)
  if(sum(colnames(group1)!=colnames(group2))>0){
    warning("gene names are different between group1 and group2")
  }
  p=rep(1,ncol(group1))
  g1_mean=colMeans(group1)
  g2_mean=colMeans(group2)
  mean_dif=g1_mean-g2_mean
  tag=which(mean_dif!=0)
    
  if(test.method=="wilcox"){
      p[tag]=sapply(tag,function(x){wilcox.test(group1[,x], group2[,x],exact=F)$p.value})
  }
  if(test.method=="t"){
      p[tag]=sapply(tag,function(x){t.test(group1[,x], group2[,x],exact=F)$p.value})
  }
    
  #for fold change 0 sum not allowed

  if(!is.null(p_adjust_method)){p=p.adjust(p,method=p_adjust_method)}
    
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
                

#' @export                
DEG_UMI=function(group1,group2,p_threshold=NULL,log2fold_threshold=NULL,p_adjust_method='fdr',test.method=c("wilcox","t")){ #it was bonferroni
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
    
  if(test.method=="wilcox"){
      p[tag]=sapply(tag,function(x){wilcox.test(group1[,x], group2[,x],exact=F)$p.value})
  }
  if(test.method=="t"){
      p[tag]=sapply(tag,function(x){t.test(group1[,x], group2[,x],exact=F)$p.value})
  }
  
  #for fold change 0 sum not allowed
  tag=which(g1_pos | g2_pos)
  foldchange[tag]=sapply(tag,function(x){log2((g1_mean[x]+1)/(g2_mean[x]+1))})
  foldchange[which(g1_pos>g2_pos)]=Inf
  foldchange[which(g1_pos<g2_pos)]=-Inf
  
  if(!is.null(p_adjust_method)){p=p.adjust(p,method=p_adjust_method)}
    
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
                
#' @export
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

#' @export                
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
                
#' @export
dropcol <- function(df, drop) {
  df <- df [, ! names(df) %in% drop, drop = FALSE]
  return(df)
}
                
#' @export
list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
  full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
           full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}
      
      
#' @export
tiling_around_region <- function(input_bed, windowsize, windownum) {
    #generate bins around region center, for each row in bed generate 2*windownum bin
    center=round((input_bed[,2]+input_bed[,3])/2)
    temp=lapply(1:length(center),function(x){  #sapply(simplify=..)
      starts=seq(center[x]-windownum*windowsize+1,center[x]+(windownum-1)*windowsize+1,windowsize)
      ends=seq(center[x]-(windownum-1)*windowsize,center[x]+windownum*windowsize,windowsize)
      out=cbind(rep(as.character(input_bed[x,1]),windownum*2),starts,ends,rep(x,windownum*2))  
      return(out)
    })
    tiling_bed=data.frame(do.call(rbind,temp),stringsAsFactors=F)
    rm(temp)
    colnames(tiling_bed)=c("chr","sta","sto","input_id")
    tiling_bed[,2]=as.numeric(tiling_bed[,2])
    tiling_bed[,3]=as.numeric(tiling_bed[,3])
    tiling_bed[,4]=as.numeric(tiling_bed[,4])
    tiling_bed[,2][tiling_bed[,2]<=0]=0
    tiling_bed[,3][tiling_bed[,3]<=1]=1
    tiling_bed$bin_id=1:nrow(tiling_bed)
    return(tiling_bed)
}  
      
#' @export
tiling_within_region <- function(input_bed, windownum, extending_windownum, extending_windowsize) {
    #generate bins within and round regions, for each row in bed generate windownum bin within region, and extending_windownum bin extend out both side
    sta=input_bed[,2]
    sto=input_bed[,3]
    windowsize=(sto-sta)/windownum
    bin_num_per_row=windownum+extending_windownum*2
    temp=lapply(1:length(sta),function(x){  #sapply(simplify=..)
      starts=c(seq(sta[x]-extending_windownum*extending_windowsize+1,sta[x],extending_windowsize),
               seq(sta[x],sto[x]-windowsize[x],length.out = windownum),
               seq(sto[x]+1,sto[x]+(extending_windownum-1)*extending_windowsize+1,extending_windowsize))
      ends=c(seq(sta[x]-(extending_windownum-1)*extending_windowsize,sta[x],extending_windowsize),
             seq(sta[x]+windowsize[x],sto[x],length.out = windownum),
             seq(sto[x]+extending_windowsize,sto[x]+extending_windownum*extending_windowsize,extending_windowsize))
      starts=round(starts)
      ends=round(ends)
      out=cbind(rep(as.character(input_bed[x,1]),bin_num_per_row),starts,ends,rep(x,bin_num_per_row))  
      return(out)
    })
    tiling_bed=data.frame(do.call(rbind,temp),stringsAsFactors=F)
    rm(temp)
    colnames(tiling_bed)=c("chr","sta","sto","input_id")
    tiling_bed[,2]=as.numeric(tiling_bed[,2])
    tiling_bed[,3]=as.numeric(tiling_bed[,3])
    tiling_bed[,4]=as.numeric(tiling_bed[,4])
    tiling_bed[,2][tiling_bed[,2]<=0]=0
    tiling_bed[,3][tiling_bed[,3]<=1]=1
    tiling_bed$bin_id=1:nrow(tiling_bed)
    return(tiling_bed)
}

#' @export
get_tiling_matrix=function(region_bed, tiling_bed, signal_bedgraph, strand=T){
    
    #intersecting
    intersect_out=bedtools_intersect(bed1 = signal_bedgraph,bed2 = tiling_bed)
    tmp_tbl=cbind(signal_bedgraph[intersect_out$overlap[,1],4],tiling_bed$bin_id[intersect_out$overlap[,2]])
    colnames(tmp_tbl)=c("signal_intensity","bin_id")
    tmp_tbl=data.table(tmp_tbl)
    tmp_tbl=tmp_tbl[,.(signal_intensity=mean(signal_intensity)),by=bin_id]

    tiling_bed$signal_intensity=0
    tiling_bed$signal_intensity[tmp_tbl$bin_id]=tmp_tbl$signal_intensity

    tiling_matrix=matrix(tiling_bed$signal_intensity,nrow=nrow(region_bed),byrow = T)
    
    if(strand){
        negative_strand_regions=which(region_bed[,6]=="-")
        if(length(negative_strand_regions)>0){
            tiling_matrix[negative_strand_regions,]=tiling_matrix[negative_strand_regions,seq(ncol(tiling_matrix),1)]
        }
    }
    return(tiling_matrix)
}
      
SparseMatrix2Matrix=function(mat,orderByRow=T){
  #get the data frame version of the sparse matrix object mat 
  summ=Matrix::summary(mat)
  matdataframe=matrix(c(summ$i,summ$j,summ$x),nrow(summ),3)
  if(orderByRow){
    matdataframe=matdataframe[order(matdataframe[,1]),]
  }
  return(matdataframe)
}
