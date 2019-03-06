suppressMessages(require(GenomicRanges))

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

bedtools_intersect=function(bed1,bed2,overlap=T,count=F,maxgap=0,minoverlap=0,ignore_strand=T,strand_col=6){
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
                
which.colmax=function(input_mat){
    return(apply(input_mat,2,function(x){which.max(x)}))
}

matrix_Aggregate=function(mat,rowby,colby,function_to_use){
  #function_to_use: mean, median, max, min, etc
  rowby=as.numeric(factor(rowby))
  colby=as.numeric(factor(colby))
  dt=data.table(cbind(rowby,mat))
  dt=dt[,lapply(.SD, function(x){function_to_use(x,na.rm = T)}),by=rowby]
  dt=dt[order(dt$rowby),]
  rownames(dt)=dt$rowby
  dt=dt[,rowby:=NULL]
  dt=t(dt)
  dt=data.table(cbind(colby,dt))
  dt=dt[,lapply(.SD, function(x){function_to_use(x,na.rm = T)}),by=colby]
  dt=dt[order(dt$colby),]
  rownames(dt)=dt$colby
  dt=dt[,colby:=NULL]
  dt=as.matrix(t(dt))
  return(dt)
}
                
firstup <- function(x) {
    x=tolower(x)
   substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    return(x)
}

colSds=function(x){
    return(matrixStats::colSds(as.matrix(x)))
}
                
rowSds=function(x){
    return(matrixStats::rowSds(as.matrix(x)))
}
                
colMaxs=function(x){
    return(matrixStats::colMaxs(as.matrix(x)))
}
                
rowMaxs=function(x){
    return(matrixStats::rowMaxs(as.matrix(x)))
}
                
colMins=function(x){
    return(matrixStats::colMins(as.matrix(x)))
}
                
rowMins=function(x){
    return(matrixStats::rowMins(as.matrix(x)))
} 
                
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

firstup <- function(x) {
    x=tolower(x)
   substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    return(x)
}
                
which.colmax=function(input_mat){
    return(apply(input_mat,2,function(x){which.max(x)}))
}