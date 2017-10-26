#the lightTree lineage algorithm
#cellbelong is merged to cellbelongtable
#in dev
devtools::use_package('data.table')
devtools::use_package('Rtsne')
devtools::use_package('proxy')
devtools::use_package('fastcluster')
devtools::use_package('optrees')
devtools::use_package('vegan')
devtools::use_package('Matrix')

#main
LightTree_Main=function(coordinates, #the coordinate for lineage finding
                        batchidx, #library id
                        timeorder, #timepoint for each library
                        batch, #which library each cell belongs to
                        total_state_num=50, #total number of state used for calculate pesudotime
                        detectionlimit=0.1, #detection limit for whether a state is valid for a timepoint, default 10% of max state-timepoint cell number
                        blockingTime=NULL,
                        cellnum_for_small_cycle=ceiling(nrow(TagMatrix)/10),
                        small_cycle_num=50,
                        loopNumPerSmallCycle=50,
                        arrows_filter_limit=0.25
){
  #generate pesudoTime
  LightTree_PerCoordSetRet=LightTree_PerCoordSet(libTimePoint=libTimePoint, #timepoint for each sample
                                                 cellbelong=cellbelong[tagCell], #which sample each cell belongs to
                                                 coordinates=rtsne_result$Y, #the coordinate for lineage finding
                                                 total_state_num=total_state_num, #total number of state used for calculate lineage
                                                 showing_state_num=showing_state_num, #number of center showing in the plot
                                                 detectionlimit=detectionlimit, #detection limit for whether a state is valid for a timepoint, default 10% of max state-timepoint cell number
                                                 blockingTime=blockingTime,
                                                 loopNumPerSmallCycle=loopNumPerSmallCycle,
                                                 arrows_filter_limit=arrows_filter_limit)
  #MST based on PesudoTime
  MSTCT=TimeConditionedMST(coordinates,LightTree_PerCoordSetRet$pseudo_time)
}

LightTree_PerCoordSet=function(libTimePoint, #timepoint for each sample
                               cellbelong, #which sample each cell belongs to
                               coordinates, #the coordinate for lineage finding
                               total_state_num=50, #total number of state used for calculate lineage
                               blockingTime=NULL,
                               loopNumPerSmallCycle=50,
                               arrows_filter_limit=0.25,
                               displaying_gp_factor=3, #
                               tolerance_factor=3, #how many sd should the algorithm tolerate to not seperate two cluster in to two continuous group
                               detectionlimit_percent_cellnum_in_time=0.1, #detection limit for whether a state is valid for a timepoint, default 10% of cell number in given time
                               detectionlimit_ratio_to_biggest_cluster=0.2 #detection limit for whether a state is valid for a timepoint, default 20% of max state-timepoint cell number
)
  {
  ptm <- proc.time()
  pseudo_timetable=matrix(0,length(cellbelong),loopNumPerSmallCycle)
  start_end_cells=list()
  
  #create cell_dist & filtering
  celldis=LightTree_celldist_filter(coordinates,tolerance_factor)
  #time
  t=(proc.time() - ptm)[3]
  t=second_to_humanReadableTime(t)
  message(paste("celldist done. time:",t[1],"h",t[2],"m",t[3],"s"))

  #loop loopNum times
  ptm <- proc.time()
  for(i in 1:loopNumPerSmallCycle){
    temp_LightTree=LightTree_core(libTimePoint=libTimePoint,
                                  cellbelong=cellbelong,
                                  coordinates=coordinates,
                                  total_state_num=total_state_num,
                                  tolerance_factor=tolerance_factor, 
                                  detectionlimit_percent_cellnum_in_time=detectionlimit_percent_cellnum_in_time, 
                                  detectionlimit_ratio_to_biggest_cluster=detectionlimit_ratio_to_biggest_cluster,
                                  blockingTime=blockingTime,
                                  cell_dist_ret=celldis)
    start_end_cells[[i]]=temp_LightTree$timepoint_center[,c(3,6)]
    pseudo_timetable[,i]=temp_LightTree$cellbelongtable$pseudo_t
    #time
    t=(proc.time() - ptm)[3]
    t=second_to_humanReadableTime(t)
  }
  message(paste("end loop",i,"time:",t[1],"h",t[2],"m",t[3],"s"))
  pseudo_time=rowMeans(pseudo_timetable,na.rm = T)
  pseudo_timeSd=rowSds(pseudo_timetable,na.rm = T)
  
  #find the average tree
  start_end_cell_df=do.call(rbind,start_end_cells)
  #the function
  message("start tree ave")
  ave_tree_df=arrows_clustering(start_end_cell_df=start_end_cell_df,
                                coordinates=coordinates,
                                arrow_strength=NULL,
                                grouping_tolerance=displaying_gp_factor,
                                arrows_filter_limit=arrows_filter_limit)
  
  return(return_obj=list(pseudo_time=pseudo_time,
                         pseudo_timeSd=pseudo_timeSd,
                         ave_tree_df=ave_tree_df))
}



LightTree_celldist_filter=function(coordinates, #the coordinate for lineage finding, 
                                                  #a list of coordinate, or a single matrix
                                  tolerance_factor=3) #how many sd should the algorithm tolerate to not seperate two cluster in to two continuous group
{
  if(is.null(dim(coordinates))){
    if(length(coordinates)==0){
      stop("empty coordinate list")
    }
    #processing list
    cell_dis=matrix(0,nrow(coordinates[[1]]),nrow(coordinates[[1]]))
    for(i in 1:length(coordinates)){
      tmp=as.matrix(dist(coordinates[[i]]))
      scale_factor=mean(tmp)/100
      cell_dis=cell_dis+tmp/scale_factor
    }
    cell_dis=cell_dis/length(coordinates)
  }else{
    cell_dis=as.matrix(dist(coordinates))
  }

  #filter outliers based on hclust single linkage
  fit=fastcluster::hclust(as.dist(cell_dis), method='single')
  groups <- cutree(fit, h=median(fit$height)+tolerance_factor*sd(fit$height))
  #rm small groups
  gp_filter_ret=gp_filter(groups,10)
  
  return(list(gp_filter_ret=gp_filter_ret,
              cell_dis=cell_dis))
}


LightTree_core=function(libTimePoint, #timepoint for each sample
                        cellbelong, #which sample each cell belongs to
                        coordinates, #the coordinate for plotting
                        total_state_num=30, #total number of state used for plot lineage
                        blockingTime=NULL, #all states in following time can only link to state at or after blockingTime
                        tolerance_factor=3, #how many sd should the algorithm tolerate to not seperate two cluster in to two continuous group
                        displaying_gp_factor=3, 
                        detectionlimit_percent_cellnum_in_time=0.05, #detection limit for whether a state is valid for a timepoint, default 10% of cell number in given time
                        detectionlimit_ratio_to_biggest_cluster=0.2, #detection limit for whether a state is valid for a timepoint, default 20% of max state-timepoint cell number
                        cell_dist_ret #returned obj from LightTree_core_celldist
                        
){
  gp_filter_ret=cell_dist_ret$gp_filter_ret
  
  #PART I
  #generating some variable
  #generate cellbelongtable from libTimePoint and cellbelong
  cellbelongtable=data.frame(
    timepoint=libTimePoint$timeorder[match(cellbelong,libTimePoint$lib)],
    pseudo_t=rep(0,length(cellbelong)),
    timepoint_state_clustering=rep(0,length(cellbelong)),
    state_center_cell_id=rep(0,length(cellbelong)),
    lib=cellbelong
  )
  timepoint_num=max(libTimePoint$timeorder)
  
  
  #PART II
  #cell belonging
  #check global cluster num and pre-run shit with kmean
  if(is.null(total_state_num)){
    properK=round(nrow(cellbelongtable)/25)
  }else{
    properK=total_state_num
  }
  
  #Kmean
  #kmean_result = kmeans(coordinates,centers = properK,iter.max = 100)
  #clustering_result=kmean_result$cluster
  #hierachical ward
  fit=fastcluster::hclust(as.dist(cell_dist_ret$cell_dis), method='ward.D2') #clustering without removing outlier
  clustering_result <- cutree(fit, h=median(fit$height)+displaying_gp_factor*sd(fit$height))
  properK=max(clustering_result)
  cellbelongtable$timepoint_state_clustering=clustering_result
  
  cell_dist_ret$cell_dis=cell_dist_ret$cell_dis[-gp_filter_ret$rm_tag,-gp_filter_ret$rm_tag]
  filtered_cluster_result=clustering_result[-gp_filter_ret$rm_tag]
  #cluster distance
  cluster_dis_singlelink=matrix(0,properK,properK)
  for(i in 2:properK){
    for(j in 1:(i-1)){
      #calculate distance with removing outlier
      cluster_dis_singlelink[i,j]=min(cell_dist_ret$cell_dis[filtered_cluster_result==i,filtered_cluster_result==j])
    }
  }
  cluster_dis_singlelink=cluster_dis_singlelink+t(cluster_dis_singlelink)
  diag(cluster_dis_singlelink)=Inf

  #cluster clusters into super cluster, threshold determined
  MST=vegan::spantree(as.dist(cell_dist_ret$cell_dis))
  segregation_threshold=median(MST$dist)+tolerance_factor*sd(MST$dist)
  
  #find timepoint_center in each state
  #state time center
  state_time_cellnum=as.matrix(table(clustering_result,cellbelongtable$timepoint))
  colnames(state_time_cellnum)=paste0("time_",1:timepoint_num)
  rownames(state_time_cellnum)=paste0("state_",1:properK)
  
  state_time_center=matrix(0,properK,timepoint_num)
  colnames(state_time_center)=paste0("time_",1:timepoint_num)
  rownames(state_time_center)=paste0("state_",1:properK)

  center_cell_for_cluster=rep(0,properK)
  for(j in 1:properK){
    tagcell=which(clustering_result==j)#same state center across all timpoint
    if(length(tagcell)>1){
      thecenter=colMeans(coordinates[tagcell,])
      cell_center_distance=apply(coordinates[tagcell,],1,function(x){sqrt(sum((x-thecenter)^2))})
      center_cell=tagcell[which.min(cell_center_distance)]        
    }else{
      if(length(tagcell)==1){
        center_cell=tagcell
      }
    }
    center_cell_for_cluster[j]=center_cell
    #update cellbelong table
    cellbelongtable$state_center_cell_id[tagcell]=center_cell
    state_time_center[j,]=center_cell
  }
  
  #time-groups filtering/assignment
  state_time_cellnum_norm=t(apply(state_time_cellnum,1,function(x){x/sum(x)}))
  for(i in 1:timepoint_num){
    fail_cell_num=which(state_time_cellnum[,i]<max(max(state_time_cellnum[,i])*detectionlimit_ratio_to_biggest_cluster,10))
    fail_cell_ratio=which(state_time_cellnum_norm[,i]<detectionlimit_percent_cellnum_in_time)
    state_time_center[union(fail_cell_num,fail_cell_ratio),i]=0
  }
  
  #move the first moment of the cluster to the earliest moment of all cluster in groups
  chained_gps=which(cluster_dis_singlelink<segregation_threshold,arr.ind = T)
  for(x in 1:timepoint_num){
    tmp=rep(0,properK)
    non_exist_state=which(state_time_center[,x]==0)
    for(y in non_exist_state){
      brother_state=chained_gps[chained_gps[,1]==y,2]
      if(length(brother_state)>1){
        for(z in 1:length(brother_state)){
          if(state_time_center[brother_state[z],x]>0){
            tmp[y]=center_cell_for_cluster[y]
            break
          }
        }
      }
    }
    state_time_center[non_exist_state,x]=tmp[non_exist_state]
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
        precusor_state[j]=prev_states[which.min(cluster_dis_singlelink[current_state[j],prev_states])]
        precusor_time[j]=which(state_time_center[precusor_state[j],]>0)
        if(length(precusor_time[j])>1){stop(paste0("more than one starting time point for state",precusor_state[j]))}
        prec[j]=state_time_center[precusor_state[j],precusor_time[j]]
        prec_state_dis[j]=min(cluster_dis_singlelink[current_state[j],prev_states])
      }
      
      #adjust prevusor in uncertainty timepoints
      #the only uncertain condition will be if the at least one of the timepoint_state do not have precusor in the same state

      #correct the lineage by moving the lineage from previous timepoint to current time point
      #if the moving can save total distance
      if(length(current_state)>1){
        #determine the number of continuous chunk in current state
        #i.e. super cluster determination
        current_dis=cluster_dis_singlelink[current_state,current_state]
        gps=get_groups(current_dis,segregation_threshold,current_state)
        for(x_gps in 1:length(gps)){
          if(length(gps[[x_gps]])>1){
            current_state_ingps=gps[[x_gps]]
            prec_state_dis_ingps=prec_state_dis[match(gps[[x_gps]],current_state)]
            #construct distance matrix between all centers
            #for current state
            current_state_dis_mat=cluster_dis_singlelink[current_state_ingps,current_state_ingps]
            current_state_dis_mat[upper.tri(current_state_dis_mat)]=0
            diag(current_state_dis_mat)=0
            #all precusor center at prev timepoint as considered as one
            current_state_dis_mat=cbind(prec_state_dis_ingps,current_state_dis_mat)
            
            #matrix to data.frame
            current_state_dis_mat=Matrix::Matrix(current_state_dis_mat,sparse=T)
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
                prec[match(gps[[x_gps]],current_state)][current_arc[2]]=starting_cell[match(gps[[x_gps]],current_state)][current_arc[1]]
                precusor_time[match(gps[[x_gps]],current_state)][current_arc[2]]=i
                precusor_state[match(gps[[x_gps]],current_state)][current_arc[2]]=current_state_ingps[current_arc[1]]
              }
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
        current_state_dis_mat=cluster_dis_singlelink[current_state,current_state]
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
  pseudo_timePerState[pseudo_timePerState==0]=NA
  cellbelongtable$pseudo_t=pseudo_timePerState[cellbelongtable$timepoint_state_clustering]
  
  ##return / end of function
  return(return_obj=list(timepoint_center=timepoint_center,
                         cellbelongtable=cellbelongtable))
}