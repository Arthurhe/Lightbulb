#the lightTree lineage algorithm
#cellbelong is merged to cellbelongtable
#in dev

devtools::use_package('Rtsne')
devtools::use_package('proxy')
devtools::use_package('matrixStats')

#main
LightTree_Main=function(TagMatrix, #the input matrix, per row cells, per column gene
                        batchidx, #library id
                        timeorder, #timepoint for each library
                        batch, #which library each cell belongs to
                        clearTime=NULL, #which timepoint is clearly depends on previous timepoint
                        total_state_num=50, #total number of state used for plot lineage
                        showing_state_num=25, #number of center showing in the plot
                        detectionlimit=0.1, #detection limit for whether a state is valid for a timepoint, default 10% of max state-timepoint cell number
                        blockingTime=NULL,
                        cellnum_for_small_cycle=ceiling(nrow(TagMatrix)/10),
                        small_cycle_num=50,
                        loopNumPerSmallCycle=50,
                        arrows_filter_limit=0.25
){
  ptm <- proc.time()
  #variable check
  if(cellnum_for_small_cycle>nrow(TagMatrix)){
    stop("cellnum_for_small_cycle must be smaller than total cell num (i.e. nrow(TagMatrix))")
  }
  
  #generating
  libTimePoint=data.frame(lib=batchidx,timeorder=timeorder)
  cellbelong=batchidx[batch]

  LightTree_PerCoordSet_returns=list()
  message (paste("LightTree scheduled to run",small_cycle_num,"cycles"))
  for(i in 1:small_cycle_num){
    #sample cells and run TSNE
    tagCell=c()
    for(j in 1:max(batch)){
      tagCell=c(tagCell,sample(which(batch==j),ceiling(cellnum_for_small_cycle/max(batch))))
    }
    tagCell=sort(sample(tagCell,cellnum_for_small_cycle))
    rtsne_result=Rtsne::Rtsne(TagMatrix[tagCell,],dims=2,max_iter = 500)
    #the arrows for each maps are
    current_ret=LightTree_PerCoordSet(libTimePoint=libTimePoint, #timepoint for each sample
                                      cellbelong=cellbelong[tagCell], #which sample each cell belongs to
                                      coordinates=rtsne_result$Y, #the coordinate for lineage finding
                                      clearTime=clearTime, #which timepoint is clearly depends on previous timepoint
                                      total_state_num=total_state_num, #total number of state used for calculate lineage
                                      showing_state_num=showing_state_num, #number of center showing in the plot
                                      detectionlimit=detectionlimit, #detection limit for whether a state is valid for a timepoint, default 10% of max state-timepoint cell number
                                      blockingTime=blockingTime,
                                      loopNumPerSmallCycle=loopNumPerSmallCycle,
                                      arrows_filter_limit=arrows_filter_limit)
    LightTree_PerCoordSet_returns[[i]]=current_ret$ave_tree_df[,c(2,4,5)]
    #adjust the cell idx to real idx
    LightTree_PerCoordSet_returns[[i]]$center=tagCell[LightTree_PerCoordSet_returns[[i]]$center]
    LightTree_PerCoordSet_returns[[i]]$prec_center[LightTree_PerCoordSet_returns[[i]]$prec_center!=0]=tagCell[LightTree_PerCoordSet_returns[[i]]$prec_center]
    message (paste("cycle",i,"finished"))
  }
  #find the average tree
  center_prec_center_pairs=do.call(rbind,LightTree_PerCoordSet_returns)
  message ("final TSNE")
  #perform full TSNE
  rtsne_result=Rtsne::Rtsne(TagMatrix,dims=2,max_iter = 1500)
  
  #the function
  ave_tree_df=arrows_clustering(start_end_cell_df=center_prec_center_pairs[,1:2],
                                arrow_strength=center_prec_center_pairs[,3],
                                showing_state_num=showing_state_num,
                                coordinates=rtsne_result$Y,
                                arrows_filter_limit=arrows_filter_limit)
  message ("done")
  t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
  message (paste("time consumed:",t[1],"hr",t[2],"min",t[3],"s"))
  return(list(coordinate=rtsne_result$Y,
              ave_tree_df=ave_tree_df,
              center_prec_center_pairs=center_prec_center_pairs))
}


LightTree_PerCoordSet=function(libTimePoint, #timepoint for each sample
                               cellbelong, #which sample each cell belongs to
                               coordinates, #the coordinate for lineage finding
                               clearTime=NULL, #which timepoint is clearly depends on previous timepoint
                               total_state_num=50, #total number of state used for calculate lineage
                               showing_state_num=25, #number of center showing in the plot
                               detectionlimit=0.1, #detection limit for whether a state is valid for a timepoint, default 10% of max state-timepoint cell number
                               blockingTime=NULL,
                               loopNumPerSmallCycle=50,
                               arrows_filter_limit=0.25)
  {
  ##main chunk
  #loop loopNum times
  pseudo_timetable=matrix(0,length(cellbelong),loopNumPerSmallCycle)
  start_end_cells=list()
  for(i in 1:loopNumPerSmallCycle){
    temp_LightTree=LightTree_core(libTimePoint=libTimePoint,
                                  cellbelong=cellbelong,
                                  clearTime=clearTime,
                                  total_state_num=total_state_num,
                                  detectionlimit=detectionlimit,
                                  coordinates=coordinates,
                                  blockingTime=blockingTime)
    start_end_cells[[i]]=temp_LightTree$timepoint_center[,c(3,6)]
    pseudo_timetable[,i]=temp_LightTree$cellbelongtable$pseudo_t
  }
  pseudo_time=rowMeans(pseudo_timetable)
  pseudo_timeSd=matrixStats::rowSds(pseudo_timetable)
  
  #find the average tree
  start_end_cell_df=do.call(rbind,start_end_cells)
  #the function
  ave_tree_df=arrows_clustering(start_end_cell_df=start_end_cell_df,
                                arrow_strength=NULL,
                                showing_state_num=showing_state_num,
                                coordinates=coordinates,
                                arrows_filter_limit=arrows_filter_limit)
  
  return(return_obj=list(pseudo_time=pseudo_time,
                         pseudo_timeSd=pseudo_timeSd,
                         ave_tree_df=ave_tree_df))
}


LightTree_core=function(libTimePoint, #timepoint for each sample
                   cellbelong, #which sample each cell belongs to
                   coordinates, #the coordinate for lineage finding
                   clearTime=NULL, #which timepoint is clearly depends on previous timepoint
                   total_state_num=30, #total number of state used for plot lineage
                   detectionlimit=0.1, #detection limit for whether a state is valid for a timepoint, default 10% of max state-timepoint cell number
                   blockingTime=NULL
                   )
  {
  #PART I
  #generating some variable
  #generate cellbelongtable from libTimePoint and cellbelong
  cellbelongtable=data.frame(
    timepoint=libTimePoint$timeorder[match(cellbelong,libTimePoint$lib)],
    pseudo_t=rep(0,nrow(coordinates)),
    timepoint_state_kmean=rep(0,nrow(coordinates)),
    state_center_cell_id=rep(0,nrow(coordinates)),
    lib=cellbelong
  )
  timepoint_num=max(libTimePoint$timeorder)
  clearTime=union(1,clearTime)
  unclearTime=setdiff(1:timepoint_num,clearTime)
  
  cell_num=nrow(coordinates)
  
  #PART II
  #cell belonging
  time_state_num=c()
  
  #check global cluster num and pre-run shit with kmean
  if(is.null(total_state_num)){
    properK=round(nrow(cellbelongtable)/50*2)
  }else{
    properK=total_state_num
  }
  
  #Kmean
  kmean_result = kmeans(coordinates,centers = properK,iter.max = 100)
  cellbelongtable$timepoint_state_kmean=kmean_result$cluster
  
  #cluster distance
  cluster_dis=as.matrix(dist(kmean_result$centers))
  diag(cluster_dis)=Inf
  

  #find timepoint_center in each state
  #state time center
  state_time_center=matrix(0,properK,timepoint_num)
  state_time_cellnum=matrix(0,properK,timepoint_num)
  colnames(state_time_center)=paste0("time_",1:timepoint_num)
  rownames(state_time_center)=paste0("state_",1:properK)
  colnames(state_time_cellnum)=paste0("time_",1:timepoint_num)
  rownames(state_time_cellnum)=paste0("state_",1:properK)
  for(j in 1:properK){
    tagcell=which(kmean_result$cluster==j)#same state center across all timpoint
    if(length(tagcell)>1){
      thecenter=colMeans(coordinates[tagcell,])
      cell_center_distance=apply(coordinates[tagcell,],1,function(x){sqrt(sum((x-thecenter)^2))})
      center_cell=tagcell[which.min(cell_center_distance)]        
    }else{
      if(length(tagcell)==1){
        center_cell=tagcell
      }
    }
    
    #update cellbelong table
    cellbelongtable$state_center_cell_id[tagcell]=center_cell
    
    for(i in 1:timepoint_num){
      tagcell=intersect(which(kmean_result$cluster==j),which(cellbelongtable$timepoint==i))
      timepointsize=sum(cellbelongtable$timepoint==i)
      state_time_cellnum[j,i]=length(tagcell)
      state_time_center[j,i]=center_cell
    }
  }
  
  #time-groups filtering/assignment
  state_time_cellnum_norm=t(apply(state_time_cellnum,1,function(x){x/sum(x)}))
  for(i in 1:timepoint_num){
    fail_cell_num=which(state_time_cellnum[,i]<max(state_time_cellnum[,i])*detectionlimit)
    fail_cell_ratio=which(state_time_cellnum[,i]<0.05)
    state_time_center[union(fail_cell_num,fail_cell_ratio),i]=0
  }
  for(j in 1:properK){
    participating_time=which(state_time_center[j,]>0)
    if(length(participating_time)>0){
      todrop=setdiff(1:ncol(state_time_cellnum),participating_time[1])
      state_time_center[j,todrop]=0
    }
  }
  
  #build tree based on existing info
  
  #PART III
  #time point center
  timepoint=c()
  timepoint_state=c()
  timepoint_center=c()
  prec_time=c()
  prec_state=c()
  prec_center=c()
  
  for(i in 1:timepoint_num){
    #what are the current timepoint_state and center cell
    current_state=which(state_time_center[,i]>0)
    #if not new state, skip     
    if(length(current_state)==0){blockingTime=setdiff(blockingTime,i);next}
    #center cell for current state
    starting_cell=state_time_center[current_state,i]
    
    if(i>1){#not first time point
      #determine the prev timepoint
      if(!is.null(blockingTime)){
        if(i<=blockingTime[1]){#determine the lastblock
          lastblock=0
        }else{
          lastblock=max(blockingTime[which(blockingTime<i)])
        }
      }else{
        lastblock=0
      }
      prev_time=lastblock:(i-1)
      
      #find precusor state center
      prev_states=which(rowSums(state_time_center[,prev_time,drop=F])>0) #which state exist
      
      #determine the precursor for each timepoint_state
      prec=c()
      precusor_state=c()
      precusor_time=c()
      prec_state_dis=c()
      curr_center_cord=coordinates[starting_cell,]
      
      for(j in 1:length(current_state)){
        precusor_state[j]=prev_states[which.min(cluster_dis[current_state[j],prev_states])]
        precusor_time[j]=which(state_time_center[precusor_state[j],]>0)
        if(length(precusor_time[j])>1){stop(paste0("more than one starting time point for state",precusor_state[j]))}
        prec[j]=state_time_center[precusor_state[j],precusor_time[j]]
        prec_state_dis[j]=min(cluster_dis[current_state[j],prev_states])
      }
      
      #adjust prevusor in uncertainty timepoints
      #the only uncertain condition will be if the at least one of the timepoint_state do not have precusor in the same state
      if(i %in% unclearTime){
        #correct the lineage by moving the lineage from previous timepoint to current time point
        #if the moving can save total distance
        if(length(current_state)>1){
          #construct distance matrix between all centers
          #for current state
          current_state_dis_mat=cluster_dis[current_state,current_state]
          current_state_dis_mat[upper.tri(current_state_dis_mat)]=0
          diag(current_state_dis_mat)=0
          #all precusor center at prev timepoint as considered as one
          current_state_dis_mat=cbind(prec_state_dis,current_state_dis_mat)
          
          #matrix to data.frame
          current_state_dis_mat=Matrix(current_state_dis_mat,sparse=T)
          current_state_dis_df=SparseMatrix2Matrix(current_state_dis_mat) #to matrix, not dataframe
          current_state_dis_df[,1]=current_state_dis_df[,1]+1
          
          #construct minimum spanning tree
          MST=optrees::getMinimumSpanningTree(1:ncol(current_state_dis_mat),current_state_dis_df,
                                              algorithm="Prim",start.node=1,show.data =F,show.graph = F)
          MST$tree.arcs[,1:2]=MST$tree.arcs[,1:2]-1
          
          #change precusor by go through the MST
          for(x in 1:nrow(MST$tree.arcs)){
            current_arc=MST$tree.arcs[x,]
            #if the originating tp_state is not precusor in prev timepoint
            if(current_arc[1]!=0){
              #change the target tp_state's precusor to cooresponding precusor in current tp
              prec[current_arc[2]]=starting_cell[current_arc[1]]
              precusor_time[current_arc[2]]=i
              precusor_state[current_arc[2]]=current_state[current_arc[1]]
            }
          }
        }
      }
    }else{
      #the first time point
      prec=rep(0,length(current_state))
      precusor_time=rep(0,length(current_state))
      precusor_state=rep(0,length(current_state))
      if(length(current_state)>1){        
        #create pseudo prec_state_dis
        #using distance to center
        #the center cluster center
        centerOf1=colMeans(kmean_result$centers[current_state,])
        topick=which.min(apply(coordinates[starting_cell,],1,function(x){sum((x-centerOf1)^2)}))
        prec_state_dis=apply(coordinates[starting_cell,],1,function(x){sum((x-coordinates[starting_cell[topick],])^2)})
        
        #construct distance matrix between all centers
        #for current state
        current_state_dis_mat=cluster_dis[current_state,current_state]
        current_state_dis_mat[upper.tri(current_state_dis_mat)]=0
        diag(current_state_dis_mat)=0
        current_state_dis_mat=cbind(prec_state_dis,current_state_dis_mat)
        
        #matrix to data.frame
        current_state_dis_mat=Matrix(current_state_dis_mat,sparse=T)
        current_state_dis_df=SparseMatrix2Matrix(current_state_dis_mat) #to matrix, not dataframe
        current_state_dis_df[,1]=current_state_dis_df[,1]+1
        
        #construct minimum spanning tree
        MST=optrees::getMinimumSpanningTree(1:ncol(current_state_dis_mat),current_state_dis_df,
                                            algorithm="Prim",start.node=1,show.data =F,show.graph = F)
        MST$tree.arcs[,1:2]=MST$tree.arcs[,1:2]-1
        
        #change precusor by go through the MST
        for(x in 1:nrow(MST$tree.arcs)){
          current_arc=MST$tree.arcs[x,]
          #if the originating tp_state is not precusor in prev timepoint
          if(current_arc[1]!=0){
            #change the target tp_state's precusor to cooresponding precusor in current tp
            prec[current_arc[2]]=starting_cell[current_arc[1]]
            precusor_time[current_arc[2]]=1
            precusor_state[current_arc[2]]=current_state[current_arc[1]]
          }
        }
      }
    }
    timepoint=c(timepoint,rep(i,length(current_state)))
    timepoint_state=c(timepoint_state,current_state)
    timepoint_center=c(timepoint_center,starting_cell)
    prec_time=c(prec_time,precusor_time)
    prec_state=c(prec_state,precusor_state)
    prec_center=c(prec_center,prec)
    if(length(precusor_state)!=length(prec)){stop("i")}
  }
  timepoint_center=data.frame(timepoint=timepoint,
                               state=timepoint_state,
                               center=timepoint_center,
                               prec_time=prec_time,
                               prec_state=prec_state,
                               prec_center=prec_center)
  
  ##the core part finished
  ##starting to calculate the timepoint information
  pseudo_timePerState=rep(0,properK)
  for(i in timepoint_center$state){
    pseudo_timePerState[i]=length(States_In_Timeline(timepoint_center,i))
  }
  cellbelongtable$pseudo_t=pseudo_timePerState[cellbelongtable$timepoint_state_kmean]
  
  ##return / end of function
  return(return_obj=list(timepoint_center=timepoint_center,
                         cellbelongtable=cellbelongtable))
}