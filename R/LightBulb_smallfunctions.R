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
                           coordinates)
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
        return(sum(arrow_strength[tag_arrow[end_states[tag_arrow]==i]]))
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