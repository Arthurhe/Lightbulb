#' @export
plotByGroup=function(Exp_Seurat,
                     Group_assignment, #must match the cells order in Exp_Seurat@dr
                     Group_type=NULL, #provide the order of Group_type
                     random_order=T,
                     embedding_type="tsne",
                     plotGroupOnly=NULL, #only plot the identified group
                     backGroundGroup=NULL, #bkg group cell as grey
                     Addtext=F,
                     main_title="",
                     target_dot_number=10000,
                     tag_cell=NULL, #plot only tag cells
                     tagcol=NULL,
                     Addlegend=T,
                     rm_small_group=F){
    
    if(sum(is.na(Group_assignment))!=0){
        warning("NAs in Group_assignment")
    }

    eval(parse(text=paste0("plotting_coordinates=Exp_Seurat@dr$",embedding_type,"@cell.embeddings[,1:2]")))
    
    #select target cell
    if(is.null(tag_cell)){
        tag_cell=1:length(Group_assignment)
    }
    plotting_coordinates=plotting_coordinates[tag_cell,]
    Group_assignment=Group_assignment[tag_cell]
    
    #rm groups that's too small
    if(rm_small_group){
        gp_size=table(Group_assignment)
        gp2rm=names(gp_size)[gp_size*min(1,target_dot_number/length(Group_assignment))<10]
        if(length(gp2rm)>0){
            to_keep=which(!Group_assignment %in% gp2rm)
            Group_assignment=Group_assignment[to_keep]
            plotting_coordinates=plotting_coordinates[to_keep,]
        }
    }
    
    Group_assignment=as.character(Group_assignment)
   
    #remove not used grouptype
    if(is.null(Group_type)){
        Group_type=unique(Group_assignment)
        Group_type=Group_type[gtools::mixedorder(Group_type)]
    }else{
        Group_type=Group_type[Group_type %in% Group_assignment]
    }
    
    #color by group id
    if(is.null(tagcol)){
        if(length(Group_type)<=9){
            tagcol=RColorBrewer::brewer.pal(max(3,length(Group_type)),"Set1")
        }else{
            tagcol=colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))
            tagcol=tagcol(length(Group_type))
        }
    }
    
    #reduce dot number to target_dot_number
    tot_dot_num=nrow(plotting_coordinates)
    if(tot_dot_num>target_dot_number){
        cell_inFg=which(!Group_assignment %in% backGroundGroup)
        cell_inBkg=which(Group_assignment %in% backGroundGroup)
        
        target_percentage_Fg=2*length(cell_inFg)/(length(Group_assignment)+length(cell_inFg))
        target_percentage_bkg=length(cell_inBkg)/(length(Group_assignment)+length(cell_inFg))
        
        if(target_percentage_bkg>0.5){
            get1w_Fg = Subsample_by_group(Group_assignment[cell_inFg], round(target_dot_number*target_percentage_Fg))
            get1w_Bkg = Subsample_by_group(Group_assignment[cell_inBkg], round(target_dot_number*target_percentage_bkg))
            get1w = c(cell_inFg[get1w_Fg],cell_inBkg[get1w_Bkg])
        }else{
            get1w=Subsample_by_group(Group_assignment,target_dot_number)
        }
        plotting_coordinates=plotting_coordinates[get1w,]
        cluster_assign=Group_assignment[get1w]
    }else{
        cluster_assign=Group_assignment
    }    
    
    if(random_order){
        random_order=sample(1:length(cluster_assign),length(cluster_assign))
    }else{
        random_order=1:length(cluster_assign)
    }
    
    #calculate cluster center
    if(Addtext){
        cluster_center=matrix(0,length(Group_type),2)
        not_present=c()
        for(i in 1:length(Group_type)){
            #check exist or not
            if(sum(cluster_assign==Group_type[i],na.rm=T)==0){
                not_present=c(not_present,i)
                next
            }
            
            #check if only 1 dot
            if(sum(cluster_assign==Group_type[i],na.rm=T)<3){
                cluster_center[i,]=plotting_coordinates[which(cluster_assign==Group_type[i])[1],]
            }else{
                tagcoord=plotting_coordinates[cluster_assign==Group_type[i],]
                quick_kmean=kmeans(tagcoord,centers=2)
                center=quick_kmean$centers[which.max(table(quick_kmean$cluster)),]
                cluster_center[i,]=tagcoord[which.min(proxy::dist(tagcoord,t(center))),]            
            }
        }
    }
    
    maxx=max(plotting_coordinates[,1])
    minx=min(plotting_coordinates[,1])
    maxy=max(plotting_coordinates[,2])
    miny=min(plotting_coordinates[,2])  
    if(is.null(plotGroupOnly)){
        if(!is.null(backGroundGroup)){
            background_cell=which(cluster_assign %in% backGroundGroup)
            random_order_b=random_order[random_order %in% background_cell]
            if(length(background_cell)==0){
                warning("backGroundGroup not found in assignment")
                #no backgroud
                plot(plotting_coordinates[random_order,],pch=19,col=tagcol[match(cluster_assign,Group_type)][random_order],cex=0.5,
                     main=paste(main_title,embedding_type),ylim=c(miny,maxy),xlim=c(minx,maxx+0.3*(maxx-minx)))   
            }else{
                #plot if yes background exist
                plot(plotting_coordinates[random_order_b,],pch=19,col="gray",cex=0.5,
                     main=paste(main_title,embedding_type),ylim=c(miny,maxy),xlim=c(minx,maxx+0.3*(maxx-minx)))
                foreground_cell=which(!cluster_assign %in% backGroundGroup)
                random_order=random_order[random_order %in% foreground_cell]
                points(plotting_coordinates[random_order,],pch=19,col=tagcol[match(cluster_assign,Group_type)][random_order],cex=0.5)
                tagcol[Group_type %in% backGroundGroup]="gray"
            }
        }else{
            #no backgroud
            plot(plotting_coordinates[random_order,],pch=19,col=tagcol[match(cluster_assign,Group_type)][random_order],cex=0.5,
                 main=paste(main_title,embedding_type),ylim=c(miny,maxy),xlim=c(minx,maxx+0.3*(maxx-minx)))        
        }
        if(Addtext){text(cluster_center[,1],cluster_center[,2],label=Group_type,font=2,cex=1.2)}
        if(Addlegend){legend("topright",Group_type,bty = "n",lty=0,pch=19,col=tagcol)}
    }else{
        if(!plotGroupOnly %in% Group_type){
            stop("cell type to plot is missing in current data")
        }
        tag_cluster_num=match(plotGroupOnly,Group_type)
        plot(plotting_coordinates[!cluster_assign %in% plotGroupOnly,],pch=19,col="gray",cex=0.5,
             main=paste(main_title,embedding_type),ylim=c(miny,maxy),xlim=c(minx,maxx+0.3*(maxx-minx)))
        for(i in 1:length(plotGroupOnly)){
            points(plotting_coordinates[cluster_assign==plotGroupOnly[i],],pch=19,col=tagcol[which(Group_type==plotGroupOnly[i])],cex=0.5)
        }
        if(Addtext){text(cluster_center[tag_cluster_num,1],cluster_center[tag_cluster_num,2]
                         ,label=Group_type[tag_cluster_num],font=2,cex=1.2)}      
        if(Addlegend){legend("topright",Group_type[tag_cluster_num],bty = "n",lty=0,pch=19,col=tagcol[tag_cluster_num])}
    }
}

#' @export
plotByGroup_1by1=function(Exp_Seurat,
                          Group_assignment,
                          Group_type=NULL, #the group ordering
                          embedding_type="tsne",
                          plotGroupOnly=NULL, #the group to plot
                          main_title="",
                          target_dot_number=2500,
                          tagcol=NULL,
                          Addtext=T,
                          rm_small_group=F,
                          tag_cell=NULL
                         ){
    if(is.null(Group_type)){
        Group_type=unique(Group_assignment)
        Group_type=Group_type[gtools::mixedorder(Group_type)]
    }else{
        Group_type=Group_type[Group_type %in% Group_assignment]
    }
    
    if(is.null(plotGroupOnly)){
        plotGroupOnly=Group_type
    }
    
    for(i in 1:length(plotGroupOnly)){
        plotByGroup(Exp_Seurat,
                    Group_assignment=Group_assignment,
                    Group_type=Group_type,
                    random_order=F,
                    embedding_type=embedding_type,
                    plotGroupOnly=plotGroupOnly[i],
                    Addtext=Addtext,
                    main_title=paste0(main_title," ",plotGroupOnly[i]),
                    target_dot_number=target_dot_number,
                    tagcol=tagcol,
                    tag_cell=tag_cell,
                    rm_small_group = rm_small_group)
    }
}

#' @export
plotByGroupA_byB=function(Exp_Seurat,
                          Group_assignment_A,
                          Group_assignment_B,
                          Group_type=NULL, #the group ordering for A
                          plotGroupOnly=NULL, #the group in B to plot
                          embedding_type="tsne",
                          main_title="",
                          target_dot_number=2500,
                          tagcol=NULL,
                          Addtext=T,
                          rm_small_group=F,
                          tag_cell=NULL  #plot only tag cells
                         ){
    Group_assignment_A=as.character(Group_assignment_A)
    Group_assignment_B=as.character(Group_assignment_B)
    
    if(is.null(Group_type)){
        Group_type=unique(Group_assignment_A)
        Group_type=Group_type[gtools::mixedorder(Group_type)]
    }else{
        Group_type=Group_type[Group_type %in% Group_assignment_A]
    }
    
    if(is.null(plotGroupOnly)){
        Group_type_B=unique(Group_assignment_B[tag_cell])
        Group_type_B=Group_type_B[gtools::mixedorder(Group_type_B)]
        plotGroupOnly=Group_type_B
    }
    
    for(i in 1:length(plotGroupOnly)){
        New_group_assignment=Group_assignment_A
        New_group_assignment[!Group_assignment_B %in% plotGroupOnly[i]]="bkg"
        
        plotByGroup(Exp_Seurat,
                    Group_assignment=New_group_assignment,
                    Group_type=Group_type,
                    random_order=F,
                    embedding_type=embedding_type,
                    backGroundGroup="bkg",
                    Addtext=Addtext,
                    main_title=paste0(main_title," ",plotGroupOnly[i]),
                    target_dot_number=target_dot_number,
                    tagcol=tagcol,
                    tag_cell=tag_cell,
                    rm_small_group = rm_small_group)
    }
}

#' @export
plotByScore=function(Exp_Seurat,
                     Score_assignment,
                     random_order=T,
                     embedding_type="tsne",
                     main_title="",
                     colorSD=1,
                     colorCenter=mean(Score_assignment),
                     target_dot_number=10000,
                     tag_cell=NULL,
                     bkg_cell=NULL, #don't supply bkg if tag is empty
                     cols=NULL,
                     bkg_col="gray",
                     Addlegend=T) #tag_cell only subsample Exp_Seurat but not Score_assignment
{   
    Score_assignment=as.numeric(Score_assignment)
    eval(parse(text=paste0("plotting_coordinates=Exp_Seurat@dr$",embedding_type,"@cell.embeddings[,1:2]")))
    
    if(is.null(tag_cell)){
        bkg_cell=NULL
    }
    
    if(is.null(bkg_cell)){
        #reduce dot number to target_dot_number
        tot_dot_num=nrow(plotting_coordinates)
        if(tot_dot_num>target_dot_number){
            get1w=Subsample_by_group(Exp_Seurat@ident,target_dot_number)
        }else{
            get1w=1:tot_dot_num
        }
    }else{
        get1w=bkg_cell
    }
    
    if(is.null(tag_cell)){
        plotting_coordinates=plotting_coordinates[get1w,]
        Score_assignment=Score_assignment[get1w]
        foreground=1:length(get1w)
        background=1
    }else{
        allc=union(tag_cell,get1w)
        plotting_coordinates=plotting_coordinates[allc,]
        foreground=which(allc %in% tag_cell)
        background=which(!allc %in% tag_cell)
    }
    
    if(random_order){
        random_order=sample(1:length(foreground),length(foreground))
        foreground=foreground[random_order]
    }
     
    if(is.null(cols)){
        rbPal <- colorRampPalette(rev(RColorBrewer::brewer.pal(7,"RdBu")))
        lowerbound=colorCenter-colorSD*sd(Score_assignment)
        upperbound=colorCenter+colorSD*sd(Score_assignment)
        breaks=c(-Inf,seq(lowerbound, upperbound, length.out=25),Inf)
        cols=rbPal(26)[as.numeric(cut(Score_assignment,breaks = breaks, include.lowest=TRUE))]
    }
    
    maxx=max(plotting_coordinates[,1])
    minx=min(plotting_coordinates[,1])
    maxy=max(plotting_coordinates[,2])
    miny=min(plotting_coordinates[,2]) 
    #plot bkg
    plot(plotting_coordinates[background,],pch=19,cex=0.3,col=bkg_col,xlim=c(minx,maxx+0.3*(maxx-minx)),main=paste(main_title,embedding_type),ylim=c(miny,maxy))      
    points(plotting_coordinates[foreground,],pch=19,cex=0.5,col=cols[random_order])  
    if(Addlegend){legend("topright",legend=round(seq(lowerbound, upperbound, length.out=5),2),bty = "n",lty=0,pch=19,col=rbPal(5))}
}


