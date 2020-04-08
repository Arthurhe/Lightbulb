get_major_gene_pattern=function(lt_obj,exp_mat,lineage_id,lineage_cluster_vec){
    #determine lineage_id
    if(!is.numeric(lineage_id)){
        if(lineage_id %in% colnames(lt_obj$lineage)){
            lineage_id=which(colnames(lt_obj$lineage) == lineage_id)
        }else{
            stop("lineage_id must be numeric lineage id or the lineage name (as shown in the lineage data frame)")
        }
    }    
    
    tag_cells=which(lt_obj$tmp_data$cluster_assignment[[lineage_id]] %in% lineage_cluster_vec)
    tag_time_vec=lt_obj$time_points_percell[tag_cells]
    tag_time_vec=factor(tag_time_vec,levels=lt_obj$time_points)
    tag_exp_mat=exp_mat[tag_cells,]
    
    #clustering genes
    ptm <- proc.time()
    gene_cors=cor(tag_exp_mat)
    dis_cors=1-gene_cors
    geneTree = fastcluster::hclust(as.dist(dis_cors), method = "average")
    dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dis_cors,deepSplit = 2,
                                                pamRespectsDendro = F, minClusterSize = 30,verbose =0)
    #initiate eigengene calculation
    all_gene_cluster=1:max(dynamicMods)
    all_PCA_outs=list()
    eigen_mat=matrix(NA,nrow(tag_exp_mat),length(all_gene_cluster))
    variance_explained=matrix(0,10,length(all_gene_cluster))
    rownames(variance_explained)=paste0("PC",1:10)
    colnames(variance_explained)=all_gene_cluster
    for(i in all_gene_cluster){
        all_PCA_outs[[i]]=prcomp(tag_exp_mat[,dynamicMods %in% i],scale=T)
        variance_explained[,i]=all_PCA_outs[[i]]$sdev[1:10]^2 / sum(all_PCA_outs[[i]]$sdev^2)
        eigen_mat[,i]=all_PCA_outs[[i]]$x[,1]
    }
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    message(paste(length(all_gene_cluster)," possible patterns found: ",t[1],"hr",t[2],"min",t[3],"s"))
    
    lineage_data=list()
    lineage_data$lineage_id=lineage_id
    lineage_data$lineage_cluster_vec=lineage_cluster_vec
    lineage_data$time=tag_time_vec
    lineage_data$time_levels=lt_obj$time_points
    lineage_data$gene_groups=as.numeric(dynamicMods)
    lineage_data$all_PCA_outs=all_PCA_outs
    lineage_data$variance_explained=variance_explained
    lineage_data$eigen_mat=eigen_mat
    #lineage_data$exp_mat=tag_exp_mat
    return(lineage_data)
}

plot_major_gene_pattern=function(lineage_data){
    eigen_gene_id=paste0("eigen",1:ncol(lineage_data$eigen_mat))
    
    eigen_vec_wide=cbind(lineage_data$time,
                         as.data.frame(lineage_data$eigen_mat))
    colnames(eigen_vec_wide)[1]=c("time")

    eigen_vec_wide=data.table(eigen_vec_wide)
    eigen_vec_tall=melt(eigen_vec_wide, id.vars = c("time"),
                        measure.vars = setdiff(colnames(eigen_vec_wide),c("time")),
                        variable.name="eigengene",value.name="expression_level")

    g=ggplot(data=eigen_vec_tall,aes(x=time, y=expression_level, fill=eigengene)) + 
    geom_violin(alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    ggtitle("eigen gene expression") +
    labs(y = "relative expression") + 
    facet_wrap(~eigengene)
    
    return(g)
}

fit_lineage_spline=function(lineage_data,df=3,verbose=T){
    fitted_spline=list()
    for(i in 1:ncol(lineage_data$eigen_mat)){
        ptm <- proc.time()
        fitted_spline[[i]]=smooth.spline(as.numeric(lineage_data$time),lineage_data$eigen_mat[,i],df=df)
        
        t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
        if(verbose){
            message(paste("eigengene ",i," fitted: ",t[1],"hr",t[2],"min",t[3],"s"))
        }
    }
    return(fitted_spline)
}

find_time_by_eigengene_spline=function(exp_mat,lineage_data,fitted_spline){
    #get eigen_gene for input
    input_eigen=matrix(NA,nrow(exp_mat),ncol(lineage_data$eigen_mat))
    for(i in 1:ncol(lineage_data$eigen_mat)){
        tag_genes=which(lineage_data$gene_groups == i)
        a=scale(exp_mat[,tag_genes],
                              lineage_data$all_PCA_outs[[i]]$center,
                              lineage_data$all_PCA_outs[[i]]$scale)
        b=lineage_data$all_PCA_outs[[i]]$rotation
        input_eigen[,i]=(a %*% b[,1])[,1]
    }
    
    #define cost function
    cost_fun = function(time,eigenGene,fitted_spline){
        predicted_eigeneGene=rep(0,length(fitted_spline))
        for(i in 1:length(fitted_spline)){
            predicted_eigeneGene[i]=predict(fitted_spline[[i]]$fit,time)$y
        }
        cost=sum((predicted_eigeneGene-eigenGene)^2)
        return(cost)
    }
    
    ptm <- proc.time()
    predicted_time=rep(NA,nrow(exp_mat))
    prediction_error=rep(NA,nrow(exp_mat))
    for(i in 1:nrow(exp_mat)){
        #R optimize is not accurate
        #need to find the rough region first
        tmp_t=seq(0.5,length(lineage_data$time_levels)+0.5)
        tmp_cost=rep(NA,length(tmp_t))
        for(j in 1:length(tmp_t)){
            tmp_cost[j]=cost_fun(tmp_t[j],input_eigen[i,],fitted_spline)
        }
        min_t=tmp_t[which.min(tmp_cost)]
        
        #optimize in local
        optimize_out=optimize(cost_fun, interval=c(min_t-0.5,min_t+0.5),
                              eigenGene=input_eigen[i,],
                              fitted_spline=fitted_spline,
                              maximum = F,
                              tol=0.025)
        predicted_time[i]=optimize_out$minimum
        prediction_error[i]=optimize_out$objective 
    }
    t=second_to_humanReadableTime(round((proc.time() - ptm)[3],2))
    message(paste(nrow(exp_mat),"cells fitted: ",t[1],"hr",t[2],"min",t[3],"s"))
    
    return(list(predicted_time=predicted_time,
                prediction_error=prediction_error))
}

batch_lineage_refine=function(lt_obj,exp_mat,target_lineage,target_endpoint,df=5){
    if(length(target_lineage) != length(target_endpoint)){
        stop("length of target_lineage & target_endpoint must be identical")
    }
    
    lineage_list=batch_get_rough_lineage(lt_obj,target_lineage,target_endpoint)
    
    lineage_data_list=list()
    fitted_spline_list=list()
    predict_out=list()
    for(i in 1:length(target_lineage)){
        message(paste0("working on lineage ",target_lineage[i],"_",target_endpoint[i]))
        lineage_data_list[[i]]=get_major_gene_pattern(lt_obj,exp_mat=exp_mat,
                                                 lineage_id=target_lineage[i],lineage_cluster_vec=lineage_list[[i]])
        fitted_spline_list[[i]]=fit_lineage_spline(lineage_data_list[[i]],df=df,verbose = F)
        predict_out[[i]]=find_time_by_eigengene_spline(exp_mat,lineage_data_list[[i]],fitted_spline_list[[i]])
    }
    return(list(target_lineage=target_lineage,
                target_endpoint=target_endpoint,
                lineage_list=lineage_list,
                lineage_data_list=lineage_data_list,
                fitted_spline_list=fitted_spline_list,
                predict_out=predict_out))
}