#the lightTree lineage algorithm
#cellbelong is merged to cellbelongtable
#in dev

devtools::use_package('Rtsne')
devtools::use_package('proxy')
devtools::use_package('matrixStats')

#main
LightTree_Main=function(TagMatrix, #the input matrix, per row cells, per column gene
                        libTimePoint, #timepoint for each sample
                        cellbelong, #which sample each cell belongs to
                        clearTime=NULL, #which timepoint is clearly depends on previous timepoint
                        total_state_num=30, #total number of state used for plot lineage
                        detectionlimit=0.1, #detection limit for whether a state is valid for a timepoint, default 10% of max state-timepoint cell number
                        coordinates=NULL,  #the coordinate for lineage finding
                        blockingTime=NULL,
                        loopNum=50,
                        arrows_filter_limit=0.1
){}


LightTree_PerCoordSet=function(TagMatrix, #the input matrix, per row cells, per column gene
                               libTimePoint, #timepoint for each sample
                               cellbelong, #which sample each cell belongs to
                               clearTime=NULL, #which timepoint is clearly depends on previous timepoint
                               total_state_num=30, #total number of state used for plot lineage
                               detectionlimit=0.1, #detection limit for whether a state is valid for a timepoint, default 10% of max state-timepoint cell number
                               coordinates=NULL,  #the coordinate for lineage finding
                               blockingTime=NULL,
                               loopNum=50,
                               arrows_filter_limit=0.1
){
  #variable check and generation
  if(is.null(coordinates) & is.null(TagMatrix)){
    stop("the converted the manifold coordinates (TSNE result) or the raw expression matrix must be provided")
  }
  
  if(is.null(coordinates)){
    #do TSNE / dimension reduction to filter noise
    rtsne_result=Rtsne::Rtsne(TagMatrix,dims=3,max_iter = 2500)
    rtsne_coordinate=rtsne_result$Y
    #DDRTree_res=DDRTree(t(TagMatrix), dimensions = 2,ncenter=GapPick$K)
    coordinates=rtsne_coordinate #t(DDRTree_res$Z)
  }
  
  ##main chunk
  #loop loopNum times
  pseudo_timetable=matrix(0,length(cellbelong),loopNum)
  start_end_cells=list()
  for(i in 1:loopNum){
    temp_LightTree=LightTree_core(libTimePoint,cellbelong,clearTime,total_state_num,detectionlimit,coordinates,blockingTime)
    start_end_cells[[i]]=temp_LightTree$timepoint_center[,c(3,6)]
    pseudo_timetable[,i]=temp_LightTree$cellbelongtable$pseudo_t
  }
  pseudo_time=rowMeans(pseudo_timetable)
  pseudo_timeSd=matrixStats::rowSds(pseudo_timetable)
  
  #find the average tree
  #start is now
  #end is previous/last
  start_end_cell_df=do.call(rbind,start_end_cells)
  start_cell_posi=coordinates[start_end_cell_df[,1],]
  end_cell_posi=coordinates[start_end_cell_df[,2],]
  
  #clustering the end point
  startend_kmean=kmeans(rbind(start_cell_posi,end_cell_posi),centers=total_state_num, iter.max = 100,nstart = 10)
  start_states=startend_kmean$cluster[1:nrow(start_cell_posi)]
  end_states=rep(0,nrow(start_cell_posi))
  end_states[which(start_end_cell_df[,2]!=0)]=startend_kmean$cluster[(nrow(start_cell_posi)+1):length(startend_kmean$cluster)]
  #find cluster center (the cell)
  dist_tocenter=proxy::dist(startend_kmean$centers,coordinates)
  centercell=sapply(1:total_state_num,function(x){which.min(dist_tocenter[x,])})
  centercell_pTSd=pseudo_timeSd[centercell]
  centercell_pTSd=log2(centercell_pTSd/max(centercell_pTSd)*10)
  
  #build tree table
  ave_tree=list()
  for(i in 1:total_state_num){
    last_states=table(end_states[which(start_states==i)])
    last_states=last_states/sum(last_states)
    #amplify the state that's unstatble
    last_states_amp=last_states
    last_states_amp[names(last_states) != "0"]=last_states[names(last_states) != "0"] * centercell_pTSd[as.numeric(names(last_states))]
    last_states=last_states[last_states_amp>arrows_filter_limit]
    
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
  ave_tree_df=do.call(rbind,ave_tree)
  
  return(return_obj=list(pseudo_time=pseudo_time,
                         pseudo_timeSd=pseudo_timeSd,
                         ave_tree_df=ave_tree_df))
}


LightTree_core=function(libTimePoint, #timepoint for each sample
                   cellbelong, #which sample each cell belongs to
                   clearTime=NULL, #which timepoint is clearly depends on previous timepoint
                   total_state_num=30, #total number of state used for plot lineage
                   detectionlimit=0.1, #detection limit for whether a state is valid for a timepoint, default 10% of max state-timepoint cell number
                   coordinates=NULL,  #the coordinate for lineage finding
                   blockingTime=NULL
){
  #PART I
  #variable check
  if(is.null(coordinates)){
    stop("the converted the manifold coordinates (TSNE result) must be provided")
  }
  
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
  #cluster_dis[is.na(cluster_dis)]=Inf
  
  #state center cell
  state_center=rep(0,properK)
  for(i in 1:properK){
    thecenter=kmean_result$centers
    tagcell=which(kmean_result$cluster==i)
    cell_center_distance=apply(coordinates[tagcell,],1,function(x){sqrt(sum((x-thecenter)^2))})
    state_center[i]=tagcell[which.min(cell_center_distance)]
  }
  
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