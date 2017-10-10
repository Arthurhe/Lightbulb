

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

plotting_lineage=function(coord,timepoint_centers,cluster_label,lib_label){
    require(gridExtra)
    #create the master dataframe
    data_tab=data.frame(x=coord[,1],y=coord[,2],
                        library=factor(lib_label,levels = unique(lib_label)),
                        state=factor(cluster_label,levels=unique(cluster_label)))
    #plot the dots
    #plot the developmental arrow
    state_center_cood=data.frame(state=timepoint_centers$state,
                                 coord_x=coord[timepoint_centers$center,1],
                                 coord_y=coord[timepoint_centers$center,2],stringsAsFactors=F)
    starting_cell=timepoint_centers[which(timepoint_centers$prec_center==0),,drop=F]
    timepoint_centers=timepoint_centers[which(timepoint_centers$prec_center!=0),,drop=F]
    
    #convert the shit to coord
    starting_cell=data_tab[starting_cell$center,,drop=F]
    mapping_table=data.frame(cbind(coord[timepoint_centers$prec_center,],coord[timepoint_centers$center,]))
    colnames(mapping_table)=c("start_x","start_y","stop_x","stop_y")

    
    g1=ggplot(data_tab, aes(x, y,colour=state)) + 
      geom_point(size=1) + 
      geom_point(data=starting_cell,inherit.aes =F,aes(x, y),shape=4,colour="firebrick4",size=2.5,stroke = 3) +
      geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y), size = 0.8,inherit.aes =F,arrow = arrow(length = unit(0.2, "cm"))) +
      theme(legend.position="none")
    
    g3=ggplot(data_tab, aes(x, y,colour=state)) + 
      geom_point(size=1) + 
      geom_point(data=starting_cell,inherit.aes =F,aes(x, y),shape=4,colour="white",size=2.5,stroke = 3) +
      geom_text(data=state_center_cood,aes(coord_x, coord_y, label = state),inherit.aes =F) +
      theme(legend.position="none")

    g2=ggplot(data_tab, aes(x, y,colour=library)) + 
      geom_point(size=1) +
      geom_point(data=starting_cell,inherit.aes =F,aes(x, y),shape=4,colour="firebrick4",size=1.5,stroke = 1.6) +
      geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y), size = 0.4,inherit.aes =F,arrow = arrow(length = unit(0.1, "cm"))) +
      theme(legend.position="none")   
 
    return(list(state_col=g1,lib_col=g2,state_num=g3))
}

plotting_state_center=function(coord,timepoint_centers,cluster_label){
    #create the master dataframe
    data_tab=data.frame(x=coord[,1],y=coord[,2],
                        state=factor(cluster_label,levels=unique(cluster_label)))
    #plot the dots
    s1=ggplot(data_tab, aes(x, y,colour=state)) + 
      geom_point(size=1) + 
      geom_text(data=data.frame(coord[timepoint_centers$center,]),aes(X1, X2,label=timepoint_centers$state),inherit.aes=F) +
      theme(legend.position="none")
    
    return(s1)
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


