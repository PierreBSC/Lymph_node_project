### Warning : this script can only be used if the script for the analysis of total cells has been run before

##I)Detailled QC of the data

#A)Library size

#General distribution
pdf("Project_Weizmann/Revision_plot/QC_non_specific.pdf",width = 5.5,height = 5.5)
lib_size=colSums(data_final)
names(lib_size)=colnames(data_final)
par(las=1)
hist(log10(lib_size),xlab="Library size of the cells (Log10)",main="Library size distribution",
     n=100,xlim=c(1,5),xaxs='i',yaxs='i')
abline(v=log10(350),lwd=2,lty=2)

#Across sequencing batches
boxplot(log10(lib_size)~batch,outline=F,xlab="Sequencing batches",
        ylab="Library size of the cells (Log10)",cex.axis=0.7,cex.lab=1.3)
abline(h=log10(350),lwd=2,lty=2)

#Across clusters
lib_size_clusters=lib_size[colnames(heatmap_data)]
boxplot(log10(lib_size_clusters)~factor(r$clusters$PCA$community[colnames(heatmap_data)],
                                        levels = order_cluster),outline=F,names=1:10,
        ylab="Library size of the cells (Log10)",xlab="Cluster",
        main="Library size across clusters",cex.lab=1.3)


summary(Number_gene_detected[colnames(data_count)])

#Across replicates 

cell_proportion=table(replicates_count,r$clusters$PCA$community)
cell_proportion=cell_proportion[c(7,6,9,10,11,13,16:17,19,20,1:4),order_cluster[-length(order_cluster)]]
sample_names=rownames(cell_proportion)
lib_size_replicates=lib_size[replicates%in%sample_names]
selected_samples=factor(replicates[replicates%in%sample_names],sample_names)


color_ordered=col=alpha(c("grey","grey","salmon","salmon","darkred","darkred",
                          "cornflowerblue","cornflowerblue","darkblue","darkblue",
                          "palegreen","palegreen","darkgreen","darkgreen"),alpha = 0.7)

boxplot(log10(lib_size_replicates)~selected_samples,outline=F,
        ylab="Library size of the cells (Log10)",xlab="Replicates",
        main="Library size across replicates",cex.lab=1.3,cex.axis=0.7,
        names=NA,col=color_ordered)





#B)Number of genes detected

Number_gene_detected=colSums(data_final>0)
names(Number_gene_detected)=colnames(data_final)
hist(log10(Number_gene_detected),100,
     xlab="Number of genes detected (Log10)",main="Number of genes detected")

#Across sequencing batches
boxplot(log10(Number_gene_detected)~batch,outline=F,xlab="Sequencing batches",cex.axis=0.7,
        ylab="Number of detected genes (Log10)")

plot(log10(lib_size),log10(Number_gene_detected),cex=0.35,
     xlab="Library size (Log10)",ylab="Number of genes (Log10)",xlim=c(1.5,4.5),ylim=c(1,3.7),
     pch=16,col=alpha("black",0.25),main="Relation between library size and number of genes detected")
cor(log10(lib_size),log10(Number_gene_detected))
abline(lm(log10(Number_gene_detected)~log10(lib_size)),lwd=2,lty=2,col="red")
##Across cluster

Number_gene_detected_cluster=Number_gene_detected[colnames(heatmap_data)]
boxplot(log10(Number_gene_detected_cluster)~factor(r$clusters$PCA$community[colnames(heatmap_data)],
                                                   levels = order_cluster),outline=F,names=1:10,
        ylab="Number of genes (Log10)",xlab="Cluster",
        main="Number of genes detected  across clusters")

#C)Proportion of ERCC genes

ERCC_genes=rownames(data_final)
ERCC_genes=ERCC_genes[grepl(pattern = "ERCC",ERCC_genes)]
ERCC_UMI=colSums(data_final[ERCC_genes,])
ERCC_UMI=ERCC_UMI/lib_size

hist(ERCC_UMI*100,xlab="Percent of ERCC among total UMIs",n=100,xlim=c(0,100),main="",xaxs="i",yaxs="i")
boxplot(ERCC_UMI~batch,outline=F,xlab="Sequencing batches",cex.axis=0.7,
        ylab="Percent of ERCC among total UMIs")
boxplot(ERCC_UMI~condition,outline=F,xlab="Conditions",
        ylab="Percent of ERCC among total UMIs")

plot(ERCC_UMI,log10(lib_size))

ERCC_cluster=ERCC_UMI[colnames(heatmap_data)]
boxplot((ERCC_cluster)~factor(r$clusters$PCA$community[colnames(heatmap_data)],
                              levels = order_cluster),outline=F,names=1:10,
        ylab="Percent of ERCC among total UMIs",xlab="Cluster",
        main="Quality of the cells across clusters",ylim=c(0,1))

dev.off()



##II)Robustness of the analysis

adjusted_rand_index=function (x, y) 
{
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y)) 
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1))) 
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + 
                                                     a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

#A)Robustness toward parameter choice

list_clustering=c()
for (k in 20:50) {
  print(k)
  r$makeKnnGraph(k=k,type='PCA',center=T,distance='cosine')
  r$getKnnClusters(method=multilevel.community,type='PCA')
  list_clustering[[k]]=r$clusters$PCA$community
}

list_clustering=list_clustering[unlist(lapply(list_clustering,FUN = length))!=0]

adjusted_RI_matrix=matrix(0,nrow = length(list_clustering),ncol = length(list_clustering))
for (i in 1:length(list_clustering)) {
  for (j in 1:length(list_clustering)) {
    adjusted_RI_matrix[i,j]=adjusted_rand_index(list_clustering[[i]],list_clustering[[j]])
    
  }
}
pdf("Project_Weizmann/Revision_plot/Stability_non_specific.pdf",height = 8,width = 8)
pheatmap(adjusted_RI_matrix,cluster_cols = F,cluster_rows = F,breaks = seq(0,1,length.out = 100))
dev.off()

#B)Robustness toward algorithm choice

r$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine')
r$getKnnClusters(method=multilevel.community,type='PCA')
r$getKnnClusters(method=infomap.community,type='PCA',name="infomap",)
r$getKnnClusters(method=walktrap.community,type='PCA',name="walktrap")

comparison_method=c(adjusted_rand_index(r$clusters$PCA$community,r$clusters$PCA$infomap),
                    adjusted_rand_index(r$clusters$PCA$community,r$clusters$PCA$walktrap))

order_louvain=table(r$clusters$PCA$community)
order_louvain=order(order_louvain,decreasing = T)

order_infomap=table(r$clusters$PCA$infomap)
order_infomap=order(order_infomap,decreasing = T)


v=table(r$clusters$PCA$infomap,r$clusters$PCA$community)
v=v/rowSums(v)
v=v[order_infomap,order_louvain]


color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("white","red"))
  x=as.numeric(x)
  if (is.null(max_scale)) {
    max_scale=quantile(x,0.99,na.rm = T)
  }
  x_prime=ifelse(x>max_scale,max_scale,x)
  x_prime=x_prime/max_scale
  x_color=f(x_prime)/255
  x_color[!complete.cases(x_color),]=c(0,0,0)
  x_color=rgb(x_color)
  return(x_color)
}

palette_grey_red=color_convertion(seq(0,1,length.out = 100))

pdf("Project_Weizmann/Revision_plot/Stability_infomap_non_specific.pdf",height = 8,width = 8)
pheatmap(v,color = palette_grey_red,treeheight_col = 0,treeheight_row = 0,
         show_rownames = F,show_colnames = F)
dev.off()

##III)Studying the Lymphocyte clusters : what are these cells 

#A)Constructing the heatmap 

table(r$clusters$PCA$infomap[(r$clusters$PCA$community==2)])
#T ? : 3
#ILC2 : 11
View(r$diffgenes$PCA$infomap$`3`)
marker_T=r$diffgenes$PCA$infomap$`3`
marker_T=marker_T[order(marker_T$M,decreasing = T),]
marker_T=marker_T[marker_T$highest,]
marker_T=rownames(marker_T[1:10,])

marker_ILC2=r$diffgenes$PCA$infomap$`11`
marker_ILC2=marker_ILC2[order(marker_ILC2$M,decreasing = T),]
marker_ILC2=marker_ILC2[marker_ILC2$highest,]
marker_ILC2=rownames(marker_ILC2[1:10,])

ILC_heatmap=cbind(TPM_data[c(marker_ILC2,marker_T),c(which(r$clusters$PCA$infomap==11),which(r$clusters$PCA$infomap==3))])

#B)Cleaning the heatmap

ILC_heatmap=as.matrix(ILC_heatmap)
lim_values=quantile(as.numeric((ILC_heatmap)),c(0.05,0.99))

ILC_heatmap[ILC_heatmap<lim_values[1]]=lim_values[1]
ILC_heatmap[ILC_heatmap>lim_values[2]]=lim_values[2]


#C)Ploting the heatmap
png("Project_Weizmann/Revision_plot/ILC_heatmap.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(130,40,4,4))
image(rotate(ILC_heatmap[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
abline(v=sum(r$clusters$PCA$infomap==11)/ncol(ILC_heatmap),lwd=15,lty=2)
box(which = "plot",lty = 1,lwd=20)
axis(2, at=seq(0,1,length.out=nrow((ILC_heatmap[,]))), 
     labels= ( rownames(ILC_heatmap)[nrow(ILC_heatmap):1]),
     las= 2, cex.axis=7, tick = F )
dev.off()


##VI)XCL1 : important cytokine 

NK_cells=(condition_count%in%c("Control","M.smeg day1","M.smeg day2") & r$clusters$PCA$community==6)
data_Xcl1=data.frame(Xcl1=TPM_data["Xcl1",NK_cells],
                     Condition=condition_count[NK_cells])

ggplot(data_Xcl1, aes(x=Condition, y=Xcl1,fill=Condition))  +  geom_violin(trim=TRUE,scale = "width", color="black",bw=1)+ 
  scale_y_continuous(name = "Expression log2(TPM)",limits = c(0,max(data_Xcl1[,"Xcl1"]))) + scale_fill_manual(values = c("grey","salmon","darkred")) +
  theme_classic() + ggtitle("Xcl1 expression NK cells") + theme(plot.title = element_text(size=22),axis.text=element_text(size = 15),axis.title = element_text(size = 12)) 
+ scale_x_discrete(labels=c("PBS","Ms day1","Ms day2"),name="   ")  

##VII)Validating the status of resident DC2

marker_mono=c("Fcgr1","Fcgr2b","Ly6c1","Ly6c2")
pdf("Project_Weizmann/Revision_plot/Violin_plot_monocytes_DC.pdf",width = 4.5,height = 3.5)
for (k in c(marker_mono)) {
  print(violin_plot_gene(k))
}
dev.off()

##VIII)Validating the status of Neutrophils 

violin_plot_gene("Csf3r")
violin_plot_gene("Cxcr2")
violin_plot_gene("Csf3r")
violin_plot_gene("Cxcr2")


