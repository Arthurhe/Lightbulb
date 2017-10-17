#all the random small functions
#in dev

devtools::use_package('Matrix')

SparseMatrix2Dataframe=function(mat,orderByRow=T){
  #get the data frame version of the sparse matrix object mat 
  summ=Matrix::summary(mat)
  matdataframe=data.frame(i=summ$i,j=summ$j,val=summ$x)
  if(orderByRow){
    matdataframe=matdataframe[order(matdataframe$i),]
  }
  return(matdataframe)
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

States_In_Timeline=function(timepoint_centers,starting_state,stoping_state=0){
  prec_state=timepoint_centers$prec_state[which(timepoint_centers$state==starting_state)]
  timeline=starting_state
  i=2
  while(prec_state!=stoping_state){
    timeline[i]=prec_state
    prec_state=timepoint_centers$prec_state[which(timepoint_centers$state==prec_state)]
    i=i+1
  }
  return(rev(timeline))
}

arrows_clustering=function(start_end_cell_df,
                           arrow_strength=NULL,
                           showing_state_num,
                           coordinates,
                           arrows_filter_limit=0.1)
{
  #start is now
  #end is previous/last
  start_cell_posi=coordinates[start_end_cell_df[,1],]
  end_cell_posi=coordinates[start_end_cell_df[,2],]
  #clustering the end point
  startend_kmean=kmeans(rbind(start_cell_posi,end_cell_posi),centers=showing_state_num, iter.max = 100,nstart = 10)
  start_states=startend_kmean$cluster[1:nrow(start_cell_posi)]
  end_states=rep(0,nrow(start_cell_posi))
  end_states[which(start_end_cell_df[,2]!=0)]=startend_kmean$cluster[(nrow(start_cell_posi)+1):length(startend_kmean$cluster)]
  #find cluster center (the cell)
  dist_tocenter=proxy::dist(startend_kmean$centers,coordinates)
  centercell=sapply(1:showing_state_num,function(x){which.min(dist_tocenter[x,])})
  
  #build tree table
  ave_tree=list()
  for(i in 1:showing_state_num){
    if(is.null(arrow_strength)){
      #if no arrow_strength, each arrow is consider as strength 1
      last_states=table(end_states[which(start_states==i)])
    }else{
      #otherwise sum the arrow strength for particular set of interaction
      tag_arrow=which(start_states==i)
      possible_prec_state=unique(end_states[tag_arrow])
      last_states=sapply(possible_prec_state,function(x){
        return(sum(arrow_strength[tag_arrow[end_states[tag_arrow]==x]]))
      })
      names(last_states)=possible_prec_state
    }
    #get rid of the ones point toward itself
    last_states=last_states[as.numeric(names(last_states))!=i]
    #get percentage
    last_states=last_states/sum(last_states)
    #filter
    last_states=last_states[last_states>arrows_filter_limit]
    
    if(length(last_states)==0){
      warning("length = 0 for current last state, means no valid origin for this state! skipping...")
    }else{
      if(0 %in% as.numeric(names(last_states))){
        percent0=last_states["0"]
        last_states=last_states[which(names(last_states)!="0")]
        ave_tree[[i]]=data.frame(state=i,
                                 center=centercell[i],
                                 prec_state=c(0,as.numeric(names(last_states))),
                                 prec_center=c(0,centercell[as.numeric(names(last_states))]),
                                 arrow_strength=c(percent0,as.vector(last_states)))
      }else{
        ave_tree[[i]]=data.frame(state=i,
                                 center=centercell[i],
                                 prec_state=as.numeric(names(last_states)),
                                 prec_center=centercell[as.numeric(names(last_states))],
                                 arrow_strength=as.vector(last_states))
      }
    }
  }
  ave_tree_df=do.call(rbind,ave_tree)
  return(ave_tree_df)
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



#########################################################
# This program is part of the SNN-Cliq method           #
# Contact Chen Xu at UNC-Charlotte for more information.# 
#########################################################
#----- example of use------#
#data<-read.table(infile, header=TRUE, sep="\t", row.names=1);
#data<-log2(data+1)
#source('SNN.R')
#SNN(data, edge_file, k=3, distance='euclidean')
#--------------------------#


SNN<-function(coordinate, k=3, distance="euclidean"){
  # other distance options refer to dist() in R
  # data is a data.frame with row represent samples
  numSpl<-nrow(coordinate)
  x<-dist(coordinate, distance, diag=TRUE, upper=TRUE)
  x<-as.matrix(x)
  IDX<-t(apply(x,1,order)[1:k,]) # knn list
  
  edges<-list()              # SNN graph
  for (i in 1:numSpl){
    j<-i
    while (j<numSpl){
      j<-j+1
      shared<-intersect(IDX[i,], IDX[j,])
      if(length(shared)>0){			
        s<-k-0.5*(match(shared, IDX[i,])+match(shared, IDX[j,]))
        strength<-max(s)
        if (strength>0)
          edges<-rbind(edges, c(i,j,strength))
      }	
    }
  }
  write.table(edges, "SNN_edges_temp", quote=FALSE, sep='\t',col.names=FALSE,row.names=FALSE)
  system(paste0("/home/ahe/SNN_cliq.py -i SNN_edges_temp -o SNN_out_temp"))
  gps=data.table::fread("SNN_out_temp",header=F,data.table=F)
  gps=gps[,1]
  return(list(gps=gps,edges=edges))
}
    
DEG_wilcox_UMI=function(group1,group2){
  p=rep(1,ncol(group1))
  foldchange=rep(1,ncol(group1))
  g1_pos=colSums(group1)>0
  g2_pos=colSums(group2)>0
  tag=which(g1_pos & g2_pos)
  p[tag]=sapply(tag,function(x){wilcox.test(group1[,x], group2[,x],exact=F)$p.value})
  #for fold change 0 sum not allowed
  tag=which(g1_pos | g2_pos)
  foldchange[tag]=sapply(tag,function(x){log2(mean(group1[,x]))-log2(mean(group2[,x]))})
  foldchange[which(g1_pos>g2_pos)]=Inf
  foldchange[which(g1_pos<g2_pos)]=-Inf
  o=cbind(log10(p),log2(foldchange))
  colnames(o)=c("pval","foldchange")
  return(o)
}

DEG_wilcox_norm=function(group1,group2,p_threshold=NULL,z_threshold=NULL,p_adjust_method="bonferroni"){
  #p is fdr/BH corrected
  if(sum(colnames(group1)!=colnames(group2))>0){
    warning("gene names are different between group1 and group2")
  }
  p=rep(1,ncol(group1))
  delta_z=rep(1,ncol(group1))
  tag=which((colSums(group1)+colSums(group2))>0)
  p[tag]=sapply(tag,function(x){wilcox.test(group1[,x], group2[,x],exact=F)$p.value})
  p=p.adjust(p,method=p_adjust_method)
  delta_z[tag]=sapply(tag,function(x){mean(group1[,x])-mean(group2[,x])})
  #check threshold
  if(!is.null(p_threshold) & !is.null(z_threshold)){
    tag=p<p_threshold & abs(delta_z) > abs(z_threshold)
    o=data.frame(log10pval=log10(p),delta_z=delta_z,DE=tag)
  }else{
    o=data.frame(log10pval=log10(p),delta_z=delta_z)
  }
  rownames(o)=colnames(group1)

  return(o)
}

group_reassigning=function(current_gps,target_gps){
  tbl_gps=table(current_gps,target_gps)
  tbl_gps_percentage=t(apply(tbl_gps,1,function(x){round(x/sum(x)*100)}))
  new_gp_description=apply(tbl_gps_percentage,1,function(x){
    tagnum=sum(x>30)
    if(tagnum>1){
      tag=order(x,decreasing = T)[1:tagnum]
      newname=paste0(colnames(tbl_gps)[tag],"(",x[tag],"%)",collapse="_")
    }else if(tagnum==1){
      tag=order(x,decreasing = T)[1]
      newname=paste0(colnames(tbl_gps)[tag],"(",x[tag],"%)")
    }else{
      tag=order(x,decreasing = T)[1]
      newname=paste0(colnames(tbl_gps)[tag],"(",x[tag],"%)_mix")
    }
    return(newname)
  })
  new_gp_assignment=apply(tbl_gps_percentage,1,function(x){
    colnames(tbl_gps)[which.max(x)]
  })
  return(obj=list(gp_assignment=new_gp_assignment,gp_description=new_gp_description,gp_table=tbl_gps))
}
