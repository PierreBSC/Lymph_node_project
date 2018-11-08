##This script is for non commerical use only
##For any question : pierre dot bost at pasteur dot fr

##AFA : Automated Flow Analysis 
##This script is based on Flow Core packge and makes it more compact
##The index sorting data extraction script has been designed for the extraction of MARS-seq plate data


library(flowCore)
library(flowViz)
library(FNN)
library(igraph)
library(pheatmap)


AFA_pipe_load=function(loading_dir="",
                       metadata_path=NULL,
                       alternative_names=NULL,specific_file=NULL,change_name=FALSE) {
  list_files=list.files(loading_dir,full.names = T)
  
  if (!is.null(specific_file)) {
    list_files=list_files[specific_file]
  }
  
  flow_data=read.flowSet(list_files)
  cat("Reading FCS files \n")
  cat(paste(length(list_files)," FCS files have been read \n"))
  
  
  correspondance_channel=pData(parameters(flow_data[[1]]))
  new_chanel_names=ifelse(is.na(correspondance_channel$desc),yes = correspondance_channel$name,no = correspondance_channel$desc)
  
  if (change_name) {
    colnames(flow_data)=as.character(new_chanel_names)
    colnames(flow_data)=make.names(colnames(flow_data))
    
    
  }
  
  metadata=NULL
  if (!is.null(metadata_path)) {
    metadata=read.table(metadata_path,header=T,sep="\t",row.names = 1)
    metadata=AnnotatedDataFrame(metadata)
    cat("Loading metadata \n")
    
  }
  
  if (!is.null(alternative_names)) {
    rownames(metadata)=alternative_names
    sampleNames(flow_data)=alternative_names
    cat("Changing names \n")
  }
  
  if (!is.null(metadata_path)) {
    flow_data@phenoData=metadata
    cat("Merging of the data \n")
  }
  cat("Change of the name \n")
  return(flow_data)
}


AFA_pipe_transform=function(flow_data,
                            chanel_to_normalise=NULL,
                            transfrom_benchmark=1) {
  
  if (data.class(flow_data)!="flowSet") {
    stop("The dataset is not a flowSet object. Load a corect dataset",call. = F)
  } 
  
  if (is.null(chanel_to_normalise)) {
    stop("No chanel to normalize has been selected",call. = F)
  }
  
  chanel_to_normalise=colnames(flow_data)[chanel_to_normalise]
  cat("Benchmark")
  trans=estimateLogicle(flow_data[[transfrom_benchmark]],chanel_to_normalise)
  cat(" : done \n")
  cat("Transformation")
  flow_data=transform(flow_data,trans)
  cat(" : done \n")
  
  return(flow_data)
}

AFA_pipe_plot=function(flow_data,
                       output_dir,
                       chanel_to_plot=NULL,
                       print_biplot=FALSE) {
  output_dir=paste(output_dir,"/AFA_plot/",sep = "")
  dir.create(path = output_dir )
  
  if (is.null(chanel_to_plot)) {
    chanel_to_plot=1:length(colnames(flow_data))
  }
  
  pdf(paste(output_dir,"density_plot.pdf",sep = ""),width = 8,height = 12)
  cat("Printing density plots \n")
  for (k in colnames(flow_data)[chanel_to_plot]) {
    temp_formula=as.formula(paste("~",k,sep = ""))
    print(densityplot(temp_formula,flow_data,overlap=0))
    cat(paste(k,"\n"))
  }
  dev.off()
  cat("Density plots : done \n")
  
  if (print_biplot) {
    cat("Printing biplots \n")
    pdf(paste(output_dir,"Biplot.pdf",sep = ""),width = 8,height = 8)
    for (i in colnames(flow_data)[chanel_to_plot]) {
      for (j in colnames(flow_data)[chanel_to_plot]) {
        if (i!=j) {
          temp_formula=make.formula(response = i,predictors = j )
          print(xyplot(temp_formula,flow_data,smooth = F))
          cat(paste(as.character(temp_formula),"\n"))
        }
      }
    }
    cat("Biplots : done \n")
    dev.off()
    
  }
  
  
}


AFA_pipe_clustering=function(flow_data,K=30,method="phenograph",
                             chanel_to_use=NULL,
                             condition_list=NULL) {
  
  if (is.null(chanel_to_use)) {
    chanel_to_use=1:length(colnames(flow_data)) 
  }
  cat("Pre-processing ")
  ###We merged the different flowframe into one big data.frame. 
  ##The data need to be compensated and transformed in order to have a meaningfull clustering
  merged_dataset=lapply(as.list(flow_data@frames),FUN =exprs )
  dataset_size=unlist(lapply(as.list(flow_data@frames),FUN = nrow ))
  condition=rep(condition_list,dataset_size)
  v=c()
  u=c()
  for (k in names(merged_dataset)) {
    v=rbind(v,merged_dataset[[k]])
    u=c(u,rep(k,nrow(merged_dataset[[k]])))
  }
  merged_dataset=v
  merged_conditions=u
  cat(" : done \n")
  cleaned_data=merged_dataset[,chanel_to_use]
  cleaned_data=unique(cleaned_data)
  
  cat("KNN construction")
  
  phenograph_clustering=Rphenograph(cleaned_data,k = K)
  phenograph_clustering=phenograph_clustering[[2]]$membership
  
  return(phenograph_clustering)
}

AFA_pipe_low_dim_embeding=function(flow_data,method="T-SNE",
                                   chanel_to_use=NULL,
                                   legend_vector=NULL,
                                   show_plot=F) {
  if (is.null(chanel_to_use)) {
    chanel_to_use=1:length(colnames(flow_data)) 
  }
  
  
  
  cat("Pre-processing ")
  merged_dataset=lapply(as.list(flow_data@frames),FUN =exprs )
  dataset_size=unlist(lapply(as.list(flow_data@frames),FUN = nrow ))
  v=c()
  u=c()
  for (k in names(merged_dataset)) {
    v=rbind(v,merged_dataset[[k]])
    u=c(u,rep(k,nrow(merged_dataset[[k]])))
  }
  merged_dataset=v
  merged_conditions=u
  cleaned_data=merged_dataset[,chanel_to_use]
  cleaned_data=unique(cleaned_data)
  cat(" : done \n")
  
  
  cat("Starting of T-SNE embedding \n")
  
  if (method=="T-SNE") {
    low_dim_map=Rtsne.multicore(cleaned_data,initial_dims = ncol(cleaned_data),verbose = T,max_iter=400,num_threads = 5)
    low_dim_map=low_dim_map$Y
  }
  
  if (method=="LargeVis") {
    low_dim_map=largeVis(t(cleaned_data), dim = 2, K = 20, n_trees = 50)
    low_dim_map=t(low_dim_map$coords)
  }
  
  
  
  if (is.null(legend_vector) & show_plot) {
    plot(low_dim_map,pch=16,cex=0.5,xlab='T-SNE 1',ylab="T-SNE 2",
         xaxt="n",yaxt="n",bty='n')
  }
  
  if (!is.null(legend_vector) & show_plot) {
    plot(low_dim_map,pch=16,cex=0.5,xlab='T-SNE 1',ylab="T-SNE 2",
         xaxt="n",yaxt="n",bty='n',col=string.to.colors(legend_vector))
  }
  return(low_dim_map)
}
##

#### Mean expression of each marker for each cluster, similarity between clusters 

AFA_pipe_cluster_plot=function(flow_data,clustering,output_dir,condition_vector,chanel_to_plot=NULL) {
  output_dir=paste(output_dir,"/Cluster_plot/",sep = "")
  dir.create(path = output_dir )
  
  if (is.null(chanel_to_plot)) {
    chanel_to_plot=1:length(colnames(flow_data))
  }
  
  cat("Pre-processing ")
  merged_dataset=lapply(as.list(flow_data@frames),FUN =exprs )
  dataset_size=unlist(lapply(as.list(flow_data@frames),FUN = nrow ))
  v=c()
  u=c()
  for (k in names(merged_dataset)) {
    v=rbind(v,merged_dataset[[k]])
    u=c(u,rep(k,nrow(merged_dataset[[k]])))
  }
  merged_dataset=v
  merged_conditions=u
  cleaned_data=merged_dataset[,chanel_to_plot]
  cleaned_data=unique(cleaned_data)
  cat(" : done \n")
  
  cluster_porportion=data.frame(Proportion=as.numeric(table(clustering)/length(clustering)))
  
  
  mean_marker_intensity=aggregate(cleaned_data,by=list(clustering),FUN=mean)
  mean_marker_intensity=t(mean_marker_intensity[,-1])
  colnames(mean_marker_intensity)=1:ncol(mean_marker_intensity)
  
  
  pdf(paste(output_dir,"/Cluster_similarity.pdf",sep = ""),width = 10,height = 10)
  order_cluster=pheatmap(cor(mean_marker_intensity),annotation_col = cluster_porportion)
  dev.off()
  
  
  pdf(paste(output_dir,"/Marker_intensity_cluster.pdf",sep = ""),width = 10,height = 6)
  pheatmap(mean_marker_intensity,annotation_col = cluster_porportion,cluster_cols = order_cluster$tree_col)
  dev.off()
  
  
  pdf(paste(output_dir,"/Cluster_distribution.pdf",sep = ""),width = 10,height = 6)
  pheatmap(table(condition_vector,clustering)/rowSums(table(condition_vector,clustering)),
           cluster_cols = order_cluster$tree_col)
  dev.off()
  
  pdf(paste(output_dir,"/Cluster_porportion",sep = ""),width = 10,height = 6)
  barplot(table(clustering)[order_cluster$tree_col$order])
  dev.off()
  
  
  
}

#### More complex analysis : extraction of the index sorting metadata,  

AFA_index_sorting=function(flow_data,correspondance_table) {
  Index_sorting_FACS=c()
  
  for (k in 1:length(flow_data)) {
    u=try(getIndexSort(flow_data[[k]]),silent = T) 
    if (class(u)=="try-error") {
      cat(paste("Error with the index sorting data of the plate",sampleNames(flow_data)[k]),", skipping to the next one","\n")
      next
    }
    cat(paste("Plate",sampleNames(flow_data)[k])," was read","\n")
    cell_position=paste(LETTERS[u$XLoc+1],u$YLoc+1,sep = "")
    u=u[,-c(1,ncol(u)-1,ncol(u))] ##Removing useless informations
    u$Position=cell_position
    v=correspondance_table[correspondance_table$Plate==sampleNames(flow_data)[k],]
    
    if (nrow(v)==0) {
      cat(paste("Correspondance between the position and the Cell ID can not be found for this plate, skipping to the next plate \n"))
      next
    }
    
    rownames(v)=v$Position
    v=v[u$Position,]
    v=na.omit(v)
    u$Cell_ID=v$Cell_ID
    Index_sorting_FACS=rbind(Index_sorting_FACS,u)
  }
  return(Index_sorting_FACS)
}



