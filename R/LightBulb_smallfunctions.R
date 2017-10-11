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


