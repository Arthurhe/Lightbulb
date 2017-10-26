#all the plotting function
#in dev

devtools::use_package('ggplot2')
#devtools::use_package('RColorBrewer')

plotting_lineage_states=function(coord,timepoint_centers,cluster_label,lib_label){
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
  
  
  g1=ggplot2::ggplot(data_tab, aes(x, y,colour=state)) + 
    ggplot2::geom_point(size=1) + 
    ggplot2::geom_point(data=starting_cell,inherit.aes =F,aes(x, y),shape=4,colour="firebrick4",size=2.5,stroke = 3) +
    ggplot2::geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y), size = 0.8,inherit.aes =F,arrow = arrow(length = unit(0.2, "cm"))) +
    ggplot2::theme(legend.position="none")
  
  g3=ggplot2::ggplot(data_tab, aes(x, y,colour=state)) + 
    ggplot2::geom_point(size=1) + 
    ggplot2::geom_point(data=starting_cell,inherit.aes =F,aes(x, y),shape=4,colour="white",size=2.5,stroke = 3) +
    ggplot2::geom_text(data=state_center_cood,aes(coord_x, coord_y, label = state),inherit.aes =F) +
    ggplot2::theme(legend.position="none")
  
  g2=ggplot2::ggplot(data_tab, aes(x, y,colour=library)) + 
    ggplot2::geom_point(size=1) +
    ggplot2::geom_point(data=starting_cell,inherit.aes =F,aes(x, y),shape=4,colour="firebrick4",size=1.5,stroke = 1.6) +
    ggplot2::geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y), size = 0.4,inherit.aes =F,arrow = arrow(length = unit(0.1, "cm"))) +
    ggplot2::theme(legend.position="none")   
  
  return(list(state_col=g1,lib_col=g2,state_num=g3))
}


plotting_lineage_pseudoT=function(coord,ave_tree_df,pseudo_time,pseudo_timeSd,lib_label){
  require(gridExtra)
  rbPal <- colorRampPalette(c('#2b83ba','#abdda4','#ffffbf','#fdae61','#d7191c'))
  #create the master dataframe
  data_tab=data.frame(x=coord[,1],y=coord[,2],
                      library=factor(lib_label,levels = unique(lib_label)),
                      pseudoT=pseudo_time,
                      pseudoSd=pseudo_timeSd)
  #plot the dots
  #plot the developmental arrow
  state_center_cood=data.frame(coord_x=coord[ave_tree_df$center,1],
                               coord_y=coord[ave_tree_df$center,2],stringsAsFactors=F)
  starting_cell=ave_tree_df[which(ave_tree_df$prec_center==0),,drop=F]
  ave_tree_df=ave_tree_df[which(ave_tree_df$prec_center!=0),,drop=F]
  
  #convert the shit to coord
  starting_cell=data_tab[starting_cell$center,,drop=F]
  mapping_table=data.frame(cbind(coord[ave_tree_df$prec_center,],coord[ave_tree_df$center,],ave_tree_df$arrow_strength*1))
  colnames(mapping_table)=c("start_x","start_y","stop_x","stop_y","arrow_strength")
  
  
  g1=ggplot2::ggplot(data_tab, aes(x, y,colour=pseudoT)) + 
    ggplot2::geom_point(size=1) + 
    ggplot2::geom_point(data=starting_cell,inherit.aes =F,aes(x, y),shape=4,colour="firebrick4",size=2.5,stroke = 3) +
    ggplot2::geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y),size = ave_tree_df$arrow_strength,inherit.aes =F,arrow = arrow(length = unit(0.2, "cm"))) +
    ggplot2::theme(legend.position="none") +
    scale_colour_gradientn(colours=rbPal(25))
  
  g3=ggplot2::ggplot(data_tab, aes(x, y,colour=pseudoSd)) + 
    ggplot2::geom_point(size=1) + 
    ggplot2::geom_point(data=starting_cell,inherit.aes =F,aes(x, y),shape=4,colour="firebrick4",size=2.5,stroke = 3) +
    ggplot2::geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y),size = ave_tree_df$arrow_strength,inherit.aes =F,arrow = arrow(length = unit(0.2, "cm"))) +
    ggplot2::theme(legend.position="none") +
    scale_colour_gradientn(colours=rbPal(25))
  
  g2=ggplot2::ggplot(data_tab, aes(x, y,colour=library)) + 
    ggplot2::geom_point(size=1) +
    ggplot2::geom_point(data=starting_cell,inherit.aes =F,aes(x, y),shape=4,colour="firebrick4",size=1.5,stroke = 1.6) +
    ggplot2::geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y), size = 0.4,inherit.aes =F,arrow = arrow(length = unit(0.1, "cm"))) +
    ggplot2::theme(legend.position="none")   
  
  return(list(pseudoT=g1,lib_col=g2,pseudoSd=g3))
}


plotting_state_center=function(coord,timepoint_centers,cluster_label){
  #create the master dataframe
  data_tab=data.frame(x=coord[,1],y=coord[,2],
                      state=factor(cluster_label,levels=unique(cluster_label)))
  #plot the dots
  s1=ggplot2::ggplot(data_tab, aes(x, y,colour=state)) + 
    ggplot2::geom_point(size=1) + 
    ggplot2::geom_text(data=data.frame(coord[timepoint_centers$center,]),aes(X1, X2,label=timepoint_centers$state),inherit.aes=F) +
    ggplot2::theme(legend.position="none")
  
  return(s1)
}

plot_MST_kid=function(coord,lib_label,MST_kid){
    data_tab=data.frame(x=coord[,1],y=coord[,2],library=factor(lib_label,levels = unique(lib_label)))
    mapping_table=data.frame(cbind(coord[-1,],coord[MST_kid,]))
    colnames(mapping_table)=c("start_x","start_y","stop_x","stop_y")
    g1=ggplot2::ggplot(data_tab, aes(x, y,colour=library)) + 
        ggplot2::geom_point(size=1) + 
        ggplot2::geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y), size = 0.4,inherit.aes =F) +
        ggplot2::theme(legend.position="none")
    return(g1)
}

plot_MST_mon=function(coord,lib_label,MST_mon){
    data_tab=data.frame(x=coord[,1],y=coord[,2],library=factor(lib_label,levels = unique(lib_label)))
    origin=which(MST_mon==0)
    origin_coord=coord[origin,]
    MST_mon=MST_mon[-origin]
    mapping_table=data.frame(cbind(coord[MST_mon,],coord[-origin,]))
    colnames(mapping_table)=c("start_x","start_y","stop_x","stop_y")
    g1=ggplot2::ggplot(data_tab, aes(x, y,colour=library)) + 
        ggplot2::geom_point(size=1) + 
        #ggplot2::geom_point(data=origin_coord,inherit.aes =F,aes(x, y),shape=4,colour="firebrick4",size=2.5,stroke = 3) +
        ggplot2::geom_segment(data=mapping_table,aes(x=start_x, xend=stop_x, y=start_y, yend=stop_y), size = 0.4,inherit.aes =F) +
        ggplot2::theme(legend.position="none")
    return(g1)
}

plot_gps_ggplot=function(coord,lib_label){
    data_tab=data.frame(x=coord[,1],y=coord[,2],library=factor(lib_label,levels = unique(lib_label)))
    cell_table=data.table(x=coord[,1],y=coord[,2],lib=lib_label)
    lib_table=as.data.frame(cell_table[,lapply(.SD, mean), by = lib])

    g1=ggplot2::ggplot(data_tab, aes(x, y,colour=library)) + 
        ggplot2::geom_point(size=1) + 
        ggplot2::geom_text(data=lib_table,aes(x, y, label = lib),inherit.aes =F) +
        ggplot2::theme(legend.position="none")
    return(g1)
}

plot_gps=function(coord,lib_label){
    plot(coord,col=lib_label,cex=0.8,pch=19)
    cell_table=data.table(x=coord[,1],y=coord[,2],lib=lib_label)
    lib_table=as.data.frame(cell_table[,lapply(.SD, mean), by = lib])
    text(lib_table$x,lib_table$y,lib_table$lib)
}


plot_gradient=function(coord,val){
    rbPal <- colorRampPalette(c('#2b83ba','#abdda4','#ffffbf','#fdae61','#d7191c'))
    cols=rbPal(25)[as.numeric(cut(val,breaks = 25))]
    plot(coord,col=cols,cex=0.8,pch=19)
}

