prepare_input=function(embedding,lineage=NULL,time_points,time_points_percell,cluster_num_df=NULL,
                       min_cell_cluster=30,min_cell_sub_cluster=10,max_cluster_num_per_timepoint_lineage=10,
                       global_clusters=NULL,min_cell_percentage_by_lineageTime_in_global_cluster=0.05,
                       min_cell_percentage_by_cluster_in_global_cluster=0.2,unskippable_timepoint=NULL){
    #set lineage a single lineage if NULL
    if(is.null(lineage)){
        lineage=data.frame(OnlyLineage=rep(T,length(time_points_percell)))
    }
    
    #set the rownames and colnames for cluster_num_df if the df is provided
    if(is.null(cluster_num_df)){
        cluster_num_df=matrix(NA,length(time_points),ncol(lineage))
        cluster_num_df=data.frame(cluster_num_df)
        rownames(cluster_num_df)=time_points
        colnames(cluster_num_df)=colnames(lineage)    
    }else{
        rownames(cluster_num_df)=time_points
        colnames(cluster_num_df)=colnames(lineage)    
    }
    
    if(is.null(unskippable_timepoint)){
        #unskippable_timepoint is list of vector, that contains the unskippable time point
        unskippable_timepoint=list()
        for(i in 1:ncol(lineage)){unskippable_timepoint[[i]]=c(0,1)}
    }
    names(unskippable_timepoint)=colnames(lineage)
    
    #construct obj
    lt_obj=list(embedding=embedding, #a data.frame, each col is a dimension, each row is a cell
                lineage=lineage, #a data.frame, each col is a lineage, each row is a cell
                time_points_percell=time_points_percell,
                time_points=time_points,
                global_clusters=global_clusters,
                cluster_num_df=cluster_num_df,
                unskippable_timepoint=unskippable_timepoint,
                parameters=list(min_cell_sub_cluster=min_cell_sub_cluster,
                                min_cell_cluster=min_cell_cluster,
                                max_cluster_num_per_timepoint_lineage=max_cluster_num_per_timepoint_lineage,
                                min_cell_percentage_by_lineageTime_in_global_cluster = min_cell_percentage_by_lineageTime_in_global_cluster,
                                min_cell_percentage_by_cluster_in_global_cluster = min_cell_percentage_by_cluster_in_global_cluster),
                tmp_data=list())
    
    #generate some tmp_data
    timepoint_cell_idx=list()
    for(l in 1:ncol(lt_obj$lineage)){
        timepoint_cell_idx[[l]]=list()
        for(i in 1:length(time_points)){
            timepoint_cell_idx[[l]][[i]]=which(time_points_percell==time_points[i] & lineage[,l])
            #edit unskippable_timepoint to avoid timepoint that has no cells
            cells_num=length(timepoint_cell_idx[[l]][[i]])
            if(cells_num==0){
                lt_obj$unskippable_timepoint[[l]]=setdiff(lt_obj$unskippable_timepoint[[l]],i)
            }
        }
    }
    lt_obj$tmp_data$timepoint_cell_idx=timepoint_cell_idx
    
    return(lt_obj)
}


lightTree=function(lt_obj){
    message("filtering outlier")
    lt_obj=outlier_filtering(lt_obj)
    
    message("clustering for each lineage at each time point")
    lt_obj=cluster_per_timepoint(lt_obj)
    
    message("organizing clustering result")
    lt_obj=cluster_organizing(lt_obj)
    
    message("sub-clustering for better distance calculation")
    lt_obj=cluster_subdivision(lt_obj)
    
    message("cluster distance calcualtion")
    lt_obj=cluster_distance(lt_obj)
    
    message("find rough tree based on clustering results")
    lt_obj=tree_search(lt_obj)
    return(lt_obj)
}



outlier_filtering=function(lt_obj){
    #filter outlier for each timepoint each lineage
    if(is.null(lt_obj$global_clusters)){
        kmean_out=kmeans(lt_obj$embedding,20,nstart=20,iter.max=200)
        lt_obj$global_clusters=kmean_out$cluster
    }
    cluster_size=table(lt_obj$global_clusters)
    
    timepoint_cell_idx=lt_obj$tmp_data$timepoint_cell_idx
    cell2remove=list()
    for(l in 1:ncol(lt_obj$lineage)){
        cell2remove[[l]]=list()
        for(i in 1:length(time_points)){
            if(length(timepoint_cell_idx[[l]][[i]]>0)){
                cluster_distribution=table(lt_obj$global_clusters[timepoint_cell_idx[[l]][[i]]])
                cluster_distribution_percentage=cluster_distribution
                for(j in 1:length(cluster_distribution_percentage)){
                    cluster_distribution_percentage[j]=cluster_distribution[j] / cluster_size[which(names(cluster_size)==names(cluster_distribution)[j])]
                }
                #remove_cells in the minor clusters
                cluster2remove=which(cluster_distribution/sum(cluster_distribution) < lt_obj$parameters$min_cell_percentage_by_lineageTime_in_global_cluster & cluster_distribution_percentage < lt_obj$parameters$min_cell_percentage_by_cluster_in_global_cluster)
                cluster2remove=names(cluster_distribution)[cluster2remove]
                cell2remove[[l]][[i]]=which(as.character(lt_obj$global_clusters[timepoint_cell_idx[[l]][[i]]]) %in% cluster2remove)

                #calculate closest neighbor distance
                cell2cell_dist=dist(lt_obj$embedding[timepoint_cell_idx[[l]][[i]],])
                cloest_neighbor_dist=rowMins(cell2cell_dist)

                #remove cells that are too far from their closest neighbor
                outlier=which((cloest_neighbor_dist-median(cloest_neighbor_dist))/mad(cloest_neighbor_dist)>3)
                cell2remove[[l]][[i]]=c(cell2remove[[l]][[i]],outlier)

                #remove outlier from list
                if(length(cell2remove[[l]][[i]]>0)){
                    timepoint_cell_idx[[l]][[i]]=timepoint_cell_idx[[l]][[i]][-cell2remove[[l]][[i]]]
                }
            }
        }
    }
    lt_obj$tmp_data$timepoint_cell_idx=timepoint_cell_idx
    lt_obj$tmp_data$cell2remove=cell2remove
    
    return(lt_obj)
}


cluster_per_timepoint=function(lt_obj){
    clustering_out=list()
    for(l in 1:ncol(lt_obj$lineage)){
        clustering_out[[l]]=list()
        ptm <- proc.time()
        for(i in 1:length(lt_obj$time_points)){
            if(length(lt_obj$tmp_data$timepoint_cell_idx[[l]][[i]])==0){
                #if there is no cell for the given lineage in given timepoint
                #set k=0 and skip
                clustering_out[[l]][[i]]=list()
                clustering_out[[l]][[i]]$k=0
            }else{
                if(is.na(lt_obj$cluster_num_df[i,l])){
                    max_cluster_num=min(round(length(lt_obj$tmp_data$timepoint_cell_idx[[l]][[i]])/lt_obj$parameters$min_cell_cluster), lt_obj$parameters$max_cluster_num_per_timepoint_lineage)
                    if(max_cluster_num>1){
                        gap_out=cluster::clusGap(lt_obj$embedding[lt_obj$tmp_data$timepoint_cell_idx[[l]][[i]],],
                                                 kmeans,nstart = 20,
                                                 iter.max = 100,
                                                 K.max=max_cluster_num,
                                                 spaceH0 ="original")
                        k=cluster::maxSE(gap_out$Tab[,3],gap_out$Tab[,4])
                    }else{
                        k=1
                    }
                }else{
                    k=min(lt_obj$cluster_num_df[i,l],round(length(lt_obj$tmp_data$timepoint_cell_idx[[l]][[i]])/lt_obj$parameters$min_cell_cluster))
                }
                
                #if more than 1 cluster
                if(k>1){
                    clustering_out[[l]][[i]]=kmeans(lt_obj$embedding[lt_obj$tmp_data$timepoint_cell_idx[[l]][[i]],],centers = k, nstart=20,iter.max=100)
                    #rm too small cells
                    tooSmall=which(clustering_out[[l]][[i]]$size < lt_obj$parameters$min_cell_cluster)
                    if(length(tooSmall)>0){
                        lt_obj$tmp_data$timepoint_cell_idx[[l]][[i]]=lt_obj$tmp_data$timepoint_cell_idx[[l]][[i]][!clustering_out[[l]][[i]]$cluster %in% tooSmall]
                        clustering_out[[l]][[i]]$cluster=clustering_out[[l]][[i]]$cluster[!clustering_out[[l]][[i]]$cluster %in% tooSmall]
                        clustering_out[[l]][[i]]$centers=clustering_out[[l]][[i]]$centers[-tooSmall,,drop=F]
                        clustering_out[[l]][[i]]$withinss=clustering_out[[l]][[i]]$withinss[-tooSmall]
                        clustering_out[[l]][[i]]$size=clustering_out[[l]][[i]]$size[-tooSmall]
                        k=k-length(tooSmall)
                    }
                }else{
                    clustering_out[[l]][[i]]=kmeans(lt_obj$embedding[lt_obj$tmp_data$timepoint_cell_idx[[l]][[i]],],centers = 1, nstart=1,iter.max=5)
                }
                clustering_out[[l]][[i]]$k=k
                lt_obj$cluster_num_df[i,l]=k
            }
        }
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        message(paste("lineage ",colnames(lt_obj$lineage)[l]," done: ",t[1],"hr",t[2],"min",t[3],"s"))
    }
    lt_obj$tmp_data$clustering_out=clustering_out
    
    return(lt_obj)
}

cluster_organizing=function(lt_obj){
    #organize the clusters as 1 to n
    total_cluster_num=rep(0,ncol(lt_obj$lineage))
    cluster_timepoint=list()
    cluster_assignment=list()
    cluster_center_embedding=list()
    
    for(l in 1:ncol(lt_obj$lineage)){
        cluster_center_list=list()
        cluster_timepoint_tmp=c()
        cluster_assignment[[l]]=rep(0,nrow(lt_obj$embedding))
        current_cluster_id=0
        
        for(i in 1:length(time_points)){
            if(lt_obj$tmp_data$clustering_out[[l]][[i]]$k!=0){
                current_cluster_id=(max(current_cluster_id)+1):(max(current_cluster_id)+lt_obj$tmp_data$clustering_out[[l]][[i]]$k)
                cluster_assignment[[l]][lt_obj$tmp_data$timepoint_cell_idx[[l]][[i]]]=current_cluster_id[lt_obj$tmp_data$clustering_out[[l]][[i]]$cluster]
                cluster_timepoint_tmp[current_cluster_id]=i
                cluster_center_list[[i]]=lt_obj$tmp_data$clustering_out[[l]][[i]]$centers
            }
        }
        cluster_timepoint[[l]]=cluster_timepoint_tmp
        cluster_center_embedding[[l]]=do.call(rbind,cluster_center_list)
        rownames(cluster_center_embedding[[l]])=1:max(current_cluster_id)
        total_cluster_num[l]=nrow(cluster_center_embedding[[l]])
    }
    
    lt_obj$tmp_data$total_cluster_num=total_cluster_num
    lt_obj$tmp_data$cluster_timepoint=cluster_timepoint
    lt_obj$tmp_data$cluster_assignment=cluster_assignment
    lt_obj$tmp_data$cluster_center_embedding=cluster_center_embedding
    
    return(lt_obj)
}

cluster_subdivision=function(lt_obj){
    #divide clusters to sub cluster for better distance calculation
    sub_cluster_belonging=list()
    sub_cluster_centers_embedding=list()
    for(l in 1:ncol(lt_obj$lineage)){
        ptm <- proc.time()
        sub_cluster_centers_list=list()
        sub_cluster_belonging_tmp=c()
        current_sub_cluster_id=0
        
        for(i in 1:lt_obj$tmp_data$total_cluster_num[l]){
            tag_cells=which(lt_obj$tmp_data$cluster_assignment[[l]]==i)
            cell_coord_tmp=lt_obj$embedding[tag_cells,]
            cell_coord_tmp[,1]=cell_coord_tmp[,1]-lt_obj$tmp_data$cluster_center_embedding[[l]][i,1]
            cell_coord_tmp[,2]=cell_coord_tmp[,2]-lt_obj$tmp_data$cluster_center_embedding[[l]][i,2]
            dist2center=sqrt(rowSums(cell_coord_tmp^2))
            outlier=(dist2center-median(dist2center))/mad(dist2center)>3
            tag_cells=tag_cells[!outlier]  
            
            if(length(tag_cells)>100){
                sub_clusters=kmeans(lt_obj$embedding[tag_cells,],centers = min(5,round(length(tag_cells)/50)), nstart=3,iter.max=50)
                sub_cluster_centers_list[[i]]=sub_clusters$centers[sub_clusters$size> lt_obj$parameters$min_cell_sub_cluster,,drop=F]
                current_sub_cluster_id=(max(current_sub_cluster_id)+1):(max(current_sub_cluster_id)+nrow(sub_cluster_centers_list[[i]]))
            }else{
                sub_cluster_centers_list[[i]]=lt_obj$tmp_data$cluster_center_embedding[[l]][i,,drop=F]
                current_sub_cluster_id=max(current_sub_cluster_id)+1
            }
            sub_cluster_belonging_tmp[current_sub_cluster_id]=i
            
        }
        sub_cluster_belonging[[l]]=sub_cluster_belonging_tmp
        sub_cluster_centers_embedding[[l]]=do.call(rbind,sub_cluster_centers_list)
        rownames(sub_cluster_centers_embedding[[l]])=1:max(current_sub_cluster_id)
        
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        message(paste("lineage ",colnames(lt_obj$lineage)[l]," sub clustering done: ",t[1],"hr",t[2],"min",t[3],"s"))
    }
    lt_obj$tmp_data$sub_cluster_belonging=sub_cluster_belonging
    lt_obj$tmp_data$sub_cluster_centers_embedding=sub_cluster_centers_embedding
    return(lt_obj)
}

cluster_distance=function(lt_obj){
    dist_df=list()
    dist_mat=list()
    for(l in 1:ncol(lt_obj$lineage)){
        #calculate cluster distance
        sub_dist_mat=as.matrix(dist(lt_obj$tmp_data$sub_cluster_centers_embedding[[l]]))
        dist_mat[[l]]=matrix(0,lt_obj$tmp_data$total_cluster_num[l],lt_obj$tmp_data$total_cluster_num[l])
        
        for(i in 2:lt_obj$tmp_data$total_cluster_num[l]){
            sub_cluster_for_i=which(lt_obj$tmp_data$sub_cluster_belonging[[l]]==i)
            for(j in 1:(i-1)){
                sub_cluster_for_j=which(lt_obj$tmp_data$sub_cluster_belonging[[l]]==j)
                dist_mat[[l]][i,j]=min(sub_dist_mat[sub_cluster_for_i,sub_cluster_for_j])
            }
        }
        #dist_mat=as.matrix(dist(cluster_center_embedding))
        dist_mat[[l]][upper.tri(dist_mat[[l]])]=0
        diag(dist_mat[[l]])=0
        dist_mat[[l]]=Matrix(dist_mat[[l]],sparse=T)
        dist_df[[l]]=SparseMatrix2Matrix(dist_mat[[l]])
    }
    
    lt_obj$tmp_data$dist_mat=dist_mat
    lt_obj$tmp_data$dist_df=dist_df
    
    return(lt_obj)
}

tree_search=function(lt_obj){
    arc_tree=list()
    for(l in 1:ncol(lt_obj$lineage)){
        arc_tbl=list()
        
        #target time point
        tag_time=which(lt_obj$cluster_num_df[,l]>0)
        if(length(tag_time)<2){stop(paste0("less than 2 time point present in lineage ",l))}
        
        dist_df=lt_obj$tmp_data$dist_df[[l]]
        
        #the first occuring time point
        i=tag_time[1]
        if(lt_obj$cluster_num_df[i,l]>1){
            current_cluster=which(lt_obj$tmp_data$cluster_timepoint[[l]]==i)
            current_dist_df=dist_df[dist_df[,1] %in% current_cluster & dist_df[,2] %in% current_cluster,, drop=F]
            current_tree=optrees::getMinimumSpanningTree(current_cluster,current_dist_df,
                                                         algorithm="Prim",start.node=1,show.data =F,show.graph = F)
            arc_tbl[[i]]=current_tree$tree.arcs
        }

        #time point 2nd to n
        for(i in tag_time[-1]){
            last_unskippable_timepoint=max(lt_obj$unskippable_timepoint[[l]][lt_obj$unskippable_timepoint[[l]]<i])
            
            past_cluster=which(lt_obj$tmp_data$cluster_timepoint[[l]]<i & lt_obj$tmp_data$cluster_timepoint[[l]] >= last_unskippable_timepoint)
            current_cluster=which(lt_obj$tmp_data$cluster_timepoint[[l]]==i)
            current_dist_df=dist_df[dist_df[,1] %in% current_cluster & dist_df[,2] %in% current_cluster,, drop=F]
            current_past_dist_df=dist_df[dist_df[,1] %in% current_cluster & dist_df[,2] %in% past_cluster,, drop=F]
            current_past_dist_df=data.table(current_past_dist_df)
            current_past_dist_df=current_past_dist_df[,.(V2=V2[which.min(V3)],V3=min(V3)),by=V1]
            current_past_dist_df=as.matrix(current_past_dist_df)
            current_past_dist_df_0=current_past_dist_df
            current_past_dist_df_0[,2]=0
            current_dist_df=rbind(current_dist_df,current_past_dist_df_0)
            current_tree=optrees::getMinimumSpanningTree(c(0,current_cluster),current_dist_df,
                                                          algorithm="Prim",start.node=0,show.data =F,show.graph = F)    
            for(j in which(current_tree$tree.arcs[,1]==0)){# column one is origin
                tag_cluster=current_tree$tree.arcs[j,2] # column two is destination
                current_tree$tree.arcs[j,1]=current_past_dist_df[which(current_past_dist_df[,1]==tag_cluster),2]
            }
            arc_tbl[[i]]=current_tree$tree.arcs
        }
        arc_tree[[l]]=do.call(rbind,arc_tbl)
        colnames(arc_tree[[l]])=c("origin","destination","weight")
        
        
    }
    
    lt_obj$tmp_data$arc_tree=arc_tree
    
    lt_obj$tmp_data$endpoints=list()
    lt_obj$tmp_data$startpoints=rep(0,2)
    for(l in 1:ncol(lt_obj$lineage)){        
        #get all possible endpoint
        lt_obj$tmp_data$endpoints[[l]]=setdiff(1:lt_obj$tmp_data$total_cluster_num[l],lt_obj$tmp_data$arc_tree[[l]][,1])
        
        #all possible startpoint
        lt_obj$tmp_data$startpoints[l]=setdiff(1:lt_obj$tmp_data$total_cluster_num[l],lt_obj$tmp_data$arc_tree[[l]][,2])
    }
    
    return(lt_obj)
}

embedding_conversion=function(lightree_obj,new_embedding){
    new_lightree_obj=lightree_obj
    closest_cell=list()
    for(i in 1:ncol(new_lightree_obj$lineage)){
        center_cell_dist=proxy::dist(lightree_obj$tmp_data$cluster_center_embedding[[i]],lightree_obj$embedding)
        closest_cell[[i]]=rep(0,nrow(center_cell_dist))
        for(j in 1:nrow(center_cell_dist)){
            closest_cell[[i]][j]=which.min(center_cell_dist[j,])
            new_lightree_obj$tmp_data$cluster_center_embedding[[i]][j,]=new_embedding[closest_cell[[i]][j],]
        }
    }
    new_lightree_obj$embedding=new_embedding
    return(new_lightree_obj)
}

plotting_lineage_rough=function(lt_obj,toplot_label=NULL){
    #if null toplot_label, label all th endpoints
    if(is.null(toplot_label)){
        toplot_label=rep("",nrow(lt_obj$embedding))
        #get end points 
        for(l in 1:ncol(lt_obj$lineage)){
            for(i in lt_obj$tmp_data$endpoints[[l]]){
                toplot_label[lt_obj$tmp_data$cluster_assignment[[l]]==i]=paste0(l,"_",i)
            }
        }
    }

    mapping_table_list=list()
    for(l in 1:ncol(lt_obj$lineage)){
        cluster_center=data.frame(lt_obj$tmp_data$cluster_center_embedding[[l]])

        #convert the shit to coor
        mapping_table_list[[l]]=data.frame(cbind(cluster_center[lt_obj$tmp_data$arc_tree[[l]][,1],],
                                                 cluster_center[lt_obj$tmp_data$arc_tree[[l]][,2],]))
        colnames(mapping_table_list[[l]])=c("start_x","start_y","stop_x","stop_y")
        mapping_table_list[[l]]$lineage=colnames(lt_obj$lineage)[l]
    }
    mapping_table=do.call(rbind,mapping_table_list)        
    
    data_tab=data.table(x=lt_obj$embedding[,1],y=lt_obj$embedding[,2],
                        to_label=toplot_label)
    data_tab_mean=data_tab[,.(x=mean(x),y=mean(y)),by=to_label]
    
    segment_col=colorRampPalette(c("gray15","gray85"))(ncol(lt_obj$lineage))
    segment_cols=segment_col[match(mapping_table$lineage,colnames(lt_obj$lineage))]
    
    g1=ggplot2::ggplot(data_tab, aes(x, y,colour=to_label)) + 
        ggplot2::geom_point(size=1) + 
        ggplot2::geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y), size = 0.8,inherit.aes =F,arrow = arrow(length = unit(0.2, "cm")),colour=segment_cols) +
        ggrepel::geom_text_repel(data=data_tab_mean,aes(x, y, label = to_label),inherit.aes =F) +
        ggplot2::theme(legend.position="none")

    return(g1)
}

plotting_lineage_refined=function(lt_obj,score,lineage_id=NULL,lineage_cluster_vec=NULL,embedding=NULL){
    if(is.null(embedding)){
        embedding=lt_obj$embedding
    }
    
    data_tab=data.table(x=embedding[,1],y=embedding[,2],
                        to_label=score)
    
    midpoint=median(score)
    coltouse=RColorBrewer::brewer.pal(3,"RdYlGn")[3:1]
    g1=ggplot2::ggplot(data_tab, aes(x, y,colour=to_label)) + 
        ggplot2::geom_point(size=1) +
        scale_color_gradient2(low=coltouse[1], mid=coltouse[2] ,high=coltouse[3], midpoint=midpoint) +
        ggplot2::theme(legend.position="none")
    
    if(!is.null(lineage_id)){
        #determine lineage_id
        if(!is.numeric(lineage_id)){
            if(lineage_id %in% colnames(lt_obj$lineage)){
                lineage_id=which(colnames(lt_obj$lineage) == lineage_id)
            }else{
                stop("lineage_id must be numeric lineage id or the lineage name (as shown in the lineage data frame)")
            }
        }    
        
        #organize the arc_tree to a single table
        cluster_center=data.frame(lt_obj$tmp_data$cluster_center_embedding[[lineage_id]])
        tag_arcs=which(lt_obj$tmp_data$arc_tree[[lineage_id]][,1] %in% lineage_cluster_vec &
                       lt_obj$tmp_data$arc_tree[[lineage_id]][,2] %in% lineage_cluster_vec)
        
        mapping_table=data.frame(cbind(cluster_center[lt_obj$tmp_data$arc_tree[[lineage_id]][tag_arcs,1],],
                                       cluster_center[lt_obj$tmp_data$arc_tree[[lineage_id]][tag_arcs,2],]))
        colnames(mapping_table)=c("start_x","start_y","stop_x","stop_y")
        
        g1=g1+geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y), size = 0.8,inherit.aes =F,arrow = arrow(length = unit(0.2, "cm")),col="gray25")
    }    
    
    return(g1)
}
    

get_rough_lineage=function(lt_obj,lineage_id,end_point_cluster_id){
    #determine lineage_id
    if(!is.numeric(lineage_id)){
        if(lineage_id %in% colnames(lt_obj$lineage)){
            lineage_id=which(colnames(lt_obj$lineage) == lineage_id)
        }else{
            stop("lineage_id must be numeric lineage id or the lineage name (as shown in the lineage data frame)")
        }
    }
    
    #get lineage
    lineage_vec=end_point_cluster_id
    while(!lt_obj$tmp_data$startpoints[lineage_id] %in% lineage_vec){
        row_id=which(lt_obj$tmp_data$arc_tree[[lineage_id]][,2] == lineage_vec[1])
        ori=lt_obj$tmp_data$arc_tree[[lineage_id]][row_id,1]
        lineage_vec=c(as.numeric(ori),lineage_vec)
    }
    
    return(lineage_vec)
}

batch_get_rough_lineage=function(lt_obj,lineage_id_vec,end_point_cluster_id_vec){
    if(length(lineage_id_vec) != length(end_point_cluster_id_vec)){
        stop("lineage_id_vec and end_point_cluster_id_vec must be the same length")
    }
    
    lineage_list=list()
    for(i in 1:length(lineage_id_vec)){
        lineage_list[[i]]=get_rough_lineage(lt_obj,lineage_id_vec[i],end_point_cluster_id_vec[i])
    }
    names(lineage_list)=paste0(lineage_id_vec,"_",end_point_cluster_id_vec)
    
    return(lineage_list)
}


get_cluster_center=function(embedding,cluster_vec){
    #get the index for the center cells for each cluster
    embedding_c=data.table(x=embedding[,1],
                           y=embedding[,2],
                           cluster=cluster_vec)
    embedding_c=embedding_c[,.(x=mean(x),y=mean(y)),by=cluster]
    
    center_cell=rep(0,nrow(embedding_c))
    for(i in 1:length(embedding_c$cluster)){
        tag_cells=which(cluster_vec==embedding_c$cluster[i])
        tmp=embedding[tag_cells,]
        tmp[,1]=(tmp[,1]-embedding_c$x[i])^2
        tmp[,2]=(tmp[,2]-embedding_c$y[i])^2
        dis2center=rowSums(tmp)  
        center_cell[i]=tag_cells[which.min(dis2center)]
    }
    embedding_c$center_cell=center_cell
    return(embedding_c)           
}


individual_gene_fit=function(supercell,lt_obj,lineage_refined,df){
    #time vec for predicting
    time_min=1
    time_max=length(lt_obj$time_points)
    time=seq(time_min,time_max,0.1)
    
    f1_spline=list()
    f1_given_time=list()
    gene_predicted_per_cell=list()
    gene_residue_per_cell=list()
    for(l in 1:length(lineage_refined$target_lineage)){
        ptm <- proc.time()
        tag_cells=which(lt_obj$lineage[,l])
        f1_spline[[l]]=list()
        f1_given_time[[l]]=matrix(NA,length(time),ncol(supercell))
        gene_predicted_per_cell[[l]]=matrix(NA,length(tag_cells),ncol(supercell))
        
        for(i in 1:ncol(supercell)){
            f1_spline[[l]][[i]]=smooth.spline(lineage_refined$predict_out[[l]]$predicted_time[tag_cells],
                                              supercell[tag_cells,i],df=df)
            gene_predicted_per_cell[[l]][,i]=predict(f1_spline[[l]][[i]]$fit, lineage_refined$predict_out[[l]]$predicted_time[tag_cells])$y
            f1_given_time[[l]][,i]=predict(f1_spline[[l]][[i]]$fit,time)$y
            
        }
        names(f1_spline[[l]])=colnames(supercell)
        colnames(f1_given_time[[l]])=colnames(supercell)
        colnames(gene_predicted_per_cell[[l]])=colnames(supercell)
        
        gene_residue_per_cell[[l]]=gene_predicted_per_cell[[l]]-supercell[tag_cells,]
        
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        message(paste(ncol(supercell),"gene fitted for",colnames(lt_obj$lineage)[l],"lineage : ",t[1],"hr",t[2],"min",t[3],"s"))   
    }
    
    f1_prediction=do.call(rbind,gene_predicted_per_cell)
    f1_residue=do.call(rbind,gene_residue_per_cell)
    
    ptm <- proc.time()
    combined_time=c()
    combined_supercell=list()
    new_lineaged_label=list()
    for(l in 1:length(lineage_refined$target_lineage)){
        tag_cells=which(lt_obj$lineage[,l])
        combined_supercell[[l]]=supercell[tag_cells,]
        combined_time=c(combined_time,lineage_refined$predict_out[[l]]$predicted_time[tag_cells])
        new_lineaged_label[[l]]=rep(l,length(tag_cells))
    }
    combined_supercell=do.call(rbind,combined_supercell)
    new_lineaged_label=do.call(c,new_lineaged_label)
    
    f0_spline=list()
    f0_prediction=matrix(NA,nrow(combined_supercell),ncol(combined_supercell))
    f0_given_time=matrix(NA,length(time),ncol(combined_supercell))
    for(i in 1:ncol(combined_supercell)){
        f0_spline[[i]]=smooth.spline(combined_time,combined_supercell[,i],df=df)
        f0_prediction[,i]=predict(f0_spline[[i]]$fit, combined_time)$y
        f0_given_time[,i]=predict(f0_spline[[i]]$fit,time)$y
    }
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    message(paste(ncol(combined_supercell),"gene null hypothesis fitted:",t[1],"hr",t[2],"min",t[3],"s"))   
    
    f0_residue=f0_prediction-combined_supercell
    
    lineage_refined$gene_fit=list(fitting_time_vec=time,
                                  f1_given_time=f1_given_time,
                                  f0_given_time=f0_given_time,
                                  f1_spline=f1_spline,
                                  f1_prediction=f1_prediction,
                                  f1_residue=f1_residue,
                                  f0_spline=f0_spline,
                                  f0_prediction=f0_prediction,
                                  f0_residue=f0_residue,
                                  df=df,
                                  new_lineaged_label=new_lineaged_label,
                                  combined_time=combined_time
                                  )    
    return(lineage_refined)
}

delta_clustering=function(lineage_refined,supercell,lineage1=1,lineage2=2){
    delta_matrix=matrix(0,length(lineage_refined$gene_fit$fitting_time_vec),ncol(supercell))

    for(i in 1:ncol(supercell)){
        delta_matrix[,i]=lineage_refined$gene_fit$f1_given_time[[lineage1]][,i]-lineage_refined$gene_fit$f1_given_time[[lineage2]][,i]
    }

    #delta correlation
    delta_cor=cor(delta_matrix)

    #clustering of gene
    delta_dis=1-delta_cor
    geneTree = fastcluster::hclust(as.dist(delta_dis), method = "average")
    dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = delta_dis,deepSplit = F,
                                                pamRespectsDendro = T, minClusterSize = 30,verbose =0)
    all_gene_cluster=1:max(dynamicMods)

    all_PCA_outs=list()
    eigen_mat=matrix(NA,nrow(delta_matrix),length(all_gene_cluster))
    variance_explained=matrix(0,10,length(all_gene_cluster))
    rownames(variance_explained)=paste0("PC",1:10)
    colnames(variance_explained)=all_gene_cluster
    for(i in all_gene_cluster){
        all_PCA_outs[[i]]=prcomp(delta_matrix[,dynamicMods %in% i],scale=T)
        variance_explained[,i]=all_PCA_outs[[i]]$sdev[1:10]^2 / sum(all_PCA_outs[[i]]$sdev^2)
        eigen_mat[,i]=all_PCA_outs[[i]]$x[,1]
    }
    
    delta_sum=colSums(abs(delta_matrix))
    delta_order=order(delta_sum,decreasing=T)
    
    lineage_refined$delta_clustering=list(delta_matrix=delta_matrix,
                                          dynamicMods=dynamicMods,
                                          delta_sum=delta_sum,
                                          delta_order=delta_order,
                                          eigen_mat=eigen_mat,
                                          variance_explained=variance_explained,
                                          all_PCA_outs=all_PCA_outs
                                          )   
    return(lineage_refined)
}

plot_delta=function(lineage_refined,x_breaks=NULL,x_label=NULL){    
    eigen_vec_wide=cbind(lineage_refined$gene_fit$fitting_time_vec,
                         as.data.frame(lineage_refined$delta_clustering$eigen_mat))
    colnames(eigen_vec_wide)[1]=c("time")

    eigen_vec_wide=data.table(eigen_vec_wide)
    eigen_vec_tall=melt(eigen_vec_wide, id.vars = c("time"),
                        measure.vars = setdiff(colnames(eigen_vec_wide),c("time")),
                        variable.name="eigengene",value.name="expression_level")
    eigen_vec_tall$eigengene=paste0("cluster_",eigen_vec_tall$eigengene)
    eigen_vec_tall$eigengene=gsub("V","",eigen_vec_tall$eigengene)
    eigen_vec_tall$eigengene=factor(eigen_vec_tall$eigengene,levels=paste0("cluster_",1:ncol(lineage_refined$delta_clustering$eigen_mat)))
    
    g=ggplot(data=eigen_vec_tall,aes(x=time, y=expression_level)) + 
        geom_line() +
        geom_hline(yintercept=0, linetype="dashed", color = "red") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        ggtitle("eigen gene expression") +
        labs(y = "relative expression") + 
        facet_wrap(~eigengene)
    
    if(!is.null(x_breaks) | !is.null(x_label)){
        g= g + scale_x_continuous(breaks = x_breaks,labels = x_label)
    }
    return(g)
}
