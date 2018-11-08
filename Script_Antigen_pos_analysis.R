## Packages needed for the analysis per se
library(pagoda2)
library(Matrix)
library(uwot)


## Packages needed for the plots

library(Cairo)
library(fifer)
library(scales)
library(ggplot2)

rotate <- function(x) t(apply(x, 2, rev))

color_convertion=function(x,max_scale=NULL) {
  f <- colorRamp(c("grey","yellow","orange","red"))
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


##### I) Loading of the data

data_final=read.table("Project_Weizmann/Final_Antigen_pos/data_final.txt",row.names = 1,sep="\t")

condition=read.table("Project_Weizmann/Final_Antigen_pos/condition.txt",row.names = 1,sep="\t")
condition=as.character(condition$x)
names(condition)=colnames(data_final)

batch=read.table("Project_Weizmann/Final_Antigen_pos/batch.txt",row.names = 1,sep="\t")
batch=as.character(batch$x)
names(batch)=colnames(data_final)

replicates=read.table("Project_Weizmann/Final_Antigen_pos/raw_condition.txt",row.names = 1,sep="\t")
replicates=as.character(replicates$x)
names(replicates)=colnames(data_final)



#####II)FIltering of the data

#A)Cell filtering
lib_size=colSums(data_final)
names(lib_size)=colnames(data_final)
par(las=1,family="serif")
hist(log10(1+lib_size),xlab="Library size of the cells (log10)",main="Librarys size distribution",n=100)
abline(v=log10(500),lwd=1,col="black",lty=2)
abline(v=log10(350),lwd=2,col="red",lty=2)

##B)Gene filtering

gene_size=rowSums(data_final)
hist(log10(1+gene_size),xlab="Gene abundance (log10)",main="Gene abundance distribution",n=100,ylim=c(0,1000))
abline(v=log10(150),lwd=2,col="red",lty=2)

##C)Creation of the final dataset 

data_count=data_final[gene_size>200,lib_size>350]
data_count=as(as.matrix(data_count),"dgCMatrix")
condition_count=condition[colnames(data_count)]
replicates_count=replicates[colnames(data_count)]
batch_count=batch[colnames(data_count)]

#III)Analysis using Pagoda2 pipeline

#A)Analysis per se

r <- Pagoda2$new(data_count,log.scale=FALSE)
r$adjustVariance(plot=T,gam.k=10,verbose = T)
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3)
r$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine')

r$getKnnClusters(method=multilevel.community,type='PCA')
r$getKnnClusters(method=infomap.community,type='PCA',name="infomap")
r$getKnnClusters(method=walktrap.community,type='PCA',name="walktrap")

#B)Visualisation using UMAP : UMAP is not implemented by PAGAODA2 pipeline intitially, R implementation by jlmelville

UMAP_all_cells=umap(r$reductions$PCA,spread = 6,
                    n_neighbors = 40,verbose =T,metric = "cosine")


##Looking at clustering results projected on UMAP embedding 

plot(UMAP_all_cells,pch=16,xaxt="n",yaxt="n",bty="n",cex=0.7,
     xlab="UMAP 1",ylab="UMAP 2",main="Louvain clustering",
     col=alpha(string.to.colors(r$clusters$PCA$community),alpha = 0.4))


#IV)Heatmap generation

###A)Selection of potential cell type marker using a home made appraoch (See method section for more details)

proportion_zero=rowSums(as.matrix(data_count)==0)/ncol(data_count)
mean_expression=rowMeans(as.matrix(data_count))


### We look for genes that have an excess of zeros for their mean expression : i.e genes whose expression is concentrated in a few cells

par(las=1)
plot(log2(mean_expression),proportion_zero,pch=16,cex=0.6,col=alpha("black",0.2),
     xlab="Log2 of mean expression",ylab="Proportion of zeros",ylim=c(0,1))
proportion_zero_fit=loess.smooth(y = proportion_zero,x = log2(mean_expression),evaluation = 1000,degree = 2,span = 1/20)
lines(proportion_zero_fit,lwd=2,col="red",lty=2)

odd_score=c()
for (k in names(mean_expression)) {
  u=abs(proportion_zero_fit$x-log2(mean_expression[k]))
  corrected_x=proportion_zero_fit$x[which.min(u)]
  corrected_y=proportion_zero_fit$y[which.min(u)]
  odd_score=c(odd_score,proportion_zero[k]-(corrected_y))
}
names(odd_score)=names(mean_expression)
odd_score=(odd_score[order(odd_score,decreasing = T)])
par(las=1,family="serif")
barplot(odd_score[1:20],horiz = T,xlim=c(0,0.4),main="Top 10 outlier score")

selected_odd_genes=names(which(odd_score>(mean(odd_score)+0.5*sd(odd_score))))

#B)Selecting marker genes 

r$getDifferentialGenes(type='PCA',verbose=T,clusterType='community',upregulated.only = T,z.threshold = 3)

size_clusters=table(r$clusters$PCA$community)
order_cluster=order(size_clusters,decreasing=T)
size_clusters=(size_clusters)[order_cluster]
size_clusters=size_clusters/ncol(data_count)
order_cluster=order_cluster[size_clusters>0.05]
order_cluster=order_cluster[!order_cluster==4]

###Additional way to define the order of the cluster : ordering clusters based on hierarchical tree 
gene_mean_expression=aggregate(t(TPM_data),by=list(r$clusters$PCA$community),FUN=mean)
rownames(gene_mean_expression)=gene_mean_expression[,1]
gene_mean_expression=t(gene_mean_expression[,-1])
gene_mean_expression=gene_mean_expression[,order_cluster]

cluster_hclust=hclust(as.dist(1-cor(gene_mean_expression[r$getOdGenes(),])^2))
plot(cluster_hclust,labels=c("DCs Mig (1)","Monocytes (1)","Monocytes (2)","DCs Mig (2)",
                             "Neutrophils","Macrophages"))
order_cluster_bis=order_cluster[cluster_hclust$order]

gene_list=c()
gene_data_frame=c()

for (k in order_cluster_bis) {
  l =  (r$diffgenes$PCA$community)[[k]]
  l=l[l$highest==T,]
  if (length(grep(pattern = "ERCC",rownames(l)))!=0) {
    l=l[-grep(pattern = "ERCC",rownames(l)),]
  }
  if (length(grep(pattern = "Rp",rownames(l)))!=0) {
    l=l[-grep(pattern = "Rp",rownames(l)),]
  }
  if (length(grep(pattern = "Rik",rownames(l)))!=0) {
    l=l[-grep(pattern = "Rik",rownames(l)),]
  }
  
  if (length(grep(pattern = "Gm",rownames(l)))!=0) {
    l=l[-grep(pattern = "Gm",rownames(l)),]
  }
  l=l[order(l$M,decreasing = T),]
  l=l[intersect(rownames(l),selected_odd_genes),]
  if (nrow(l)==0) {
    next
  }
  gene_list=c(gene_list,rownames(l[1:5,]))
  gene_data_frame=rbind(gene_data_frame,data.frame(Gene=rownames(l),Score=l$Z,Cluster=k))
}

gene_list=unique(gene_list)
gene_list=intersect(gene_list,rownames(data_count))
gene_data_frame=na.omit(gene_data_frame)


#C)Computing normalised expression and concatenating gene expression into heatmap

TPM_data=t(log2(1+t(as.matrix(data_count))/lib_size[colnames(data_count)]*10^6))
TPM_data=as.matrix(TPM_data)

#Ordering clusters based on the hierarchical clustering

heatmap_data=c()
for (k in order_cluster_bis) {
  v=TPM_data[intersect(gene_list,selected_odd_genes),r$clusters$PCA$community==k]
  v=v[,sample(1:ncol(v),replace = F,size = ncol(v))]
  heatmap_data=cbind(heatmap_data,v)
  cat(paste(k,"\n"))
}

#D)Plot the heatmap

heatmap_data=as.matrix(heatmap_data)
lim_values=quantile(as.numeric((heatmap_data)),c(0.0,0.99)) ##Trimming extreme values

heatmap_data[heatmap_data<lim_values[1]]=lim_values[1]
heatmap_data[heatmap_data>lim_values[2]]=lim_values[2]

png("Project_Weizmann/Revision_plot/Heatmap_ag_pos.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(130,40,4,4))
image(rotate(heatmap_data[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(order_cluster_bis)])/sum(size_clusters[as.character(order_cluster_bis)]),lwd=10,lty=2)
axis(2, at=seq(0,1,length.out=nrow((heatmap_data[,]))), 
     labels= ( rownames(heatmap_data)[nrow(heatmap_data):1]),
     las= 2, cex.axis=7, tick = F )
dev.off()

#V)Analysis of the proportion of cells across samples

#A)Computation of the proportion

##We want to select only the most important populations

plot(1:length(size_clusters),cumsum(size_clusters),
     xlab="Number of cluster",ylab="Cumulative percentage",type="l",lwd=2)
abline(h=0.90,lwd=2,lty=2,col="red")
abline(v=7,lwd=2,lty=2,col="red")

size_clusters=table(r$clusters$PCA$community)
order_cluster=order(size_clusters,decreasing=T)
size_clusters=(size_clusters)[order_cluster]
size_clusters=size_clusters/ncol(data_count)
order_cluster=order_cluster[size_clusters>0.05]
order_cluster=order_cluster[!order_cluster==4] #ERCC cluster : low quality cluster

cell_proportion=table(replicates_count,r$clusters$PCA$community)
cell_proportion=cell_proportion[,order_cluster]

### The 2 monocytes and DC populations are mixed to improve robustness
cell_proportion=cbind(rowSums(cell_proportion[,c(1,4)]),
                      rowSums(cell_proportion[,c(2,3)]),
                      cell_proportion[,5:6])
cell_proportion=cell_proportion/rowSums(cell_proportion)*100

cell_proportion=cell_proportion[c(6:13,1,2,4,5),]

#B)Plot

black_palette=colorRampPalette(colors = c("black","white"))(9)
par(mar=c(15,2,2,2))
barplot(t(cell_proportion),las=2,col=rainbow(4))

####VI)Comparison of populations : DE analysis between the two Monocytes and DC populations

##The general approach used for DE genes is explained in details in the method section of the paper 

#A)Monocytes 

data_DE=data_count[,r$clusters$PCA$community==2 | r$clusters$PCA$community==10]
lib_size_DE=log10(colSums(data_final[,colnames(data_DE)]))
condition_DE=(r$clusters$PCA$community)[colnames(data_DE)]
batch_DE=batch_count[r$clusters$PCA$community==2 | r$clusters$PCA$community==10]

model_selection=c()
for (gene in rownames(data_DE)) {
  x=ifelse(as.numeric(data_DE[gene,])!=0,yes = 1,no = 0)
  model_logistic=glm(x~lib_size_DE+condition_DE+batch_DE,family = "binomial")
  u=summary(model_logistic)
  u=c(u$coefficients[2:3,4],u$coefficients[2:3,1])
  model_selection=rbind(model_selection,u)
}

rownames(model_selection)=rownames(data_DE)
model_selection=data.frame(model_selection)
colnames(model_selection)=c("Lib_size","Condition","Intercep","Effect")
model_selection$Condition=p.adjust(model_selection$Condition,method = "BH")
model_selection$Lib_size=p.adjust(model_selection$Lib_size,method = "BH")
model_selection$Condition=-log10(model_selection$Condition)
model_selection$Lib_size=-log10(model_selection$Lib_size)
model_selection=model_selection[order(model_selection$Condition,decreasing = T),]

mean_mono_1_expression=rowMeans(TPM_data[,r$clusters$PCA$community==2])
mean_mono_2_expression=rowMeans(TPM_data[,r$clusters$PCA$community==10])
fit_mono_logFC=loess.smooth(y = mean_mono_2_expression,x = mean_mono_1_expression,evaluation = 1000,degree = 2,span = 1/20)
corrected_mono=c()
for (k in names(mean_mono_1_expression)) {
  u=which.min(abs(fit_mono_logFC$x-(mean_mono_1_expression[k])))
  corrected_mono=c(corrected_mono,mean_mono_2_expression[k]-fit_DC_logFC$y[u])
}

Mono_valid_genes=rownames(model_selection)[model_selection$Lib_size>2] ##We only keep genes that are more likely detected if the library size is higher....

Mono_changed_genes=rownames(model_selection)[(model_selection$Condition>2 & model_selection$Lib_size>2)]
Mono_changed_genes=intersect(Mono_changed_genes,names(which(abs(corrected_mono)>2)))

plot(corrected_mono[Mono_valid_genes],model_selection[Mono_valid_genes,"Condition"],pch=16,xlim=c(-6.5,6.5),
     xlab="Monocytes (2) vs Monocytes (1) (Log2 Fold Change)",ylab="-Log10(P-value)",cex.lab=1.5,
     col=string.to.colors(Mono_valid_genes%in%Mono_changed_genes,colors = c(alpha("red",0.7),alpha("black",0.3))))

#B)DC 

data_DE=data_count[,r$clusters$PCA$community==8 | r$clusters$PCA$community==9]
lib_size_DE=log10(colSums(data_final[,colnames(data_DE)]))
condition_DE=(r$clusters$PCA$community)[colnames(data_DE)]
batch_DE=batch_count[r$clusters$PCA$community==8 | r$clusters$PCA$community==9]

model_selection=c()
for (gene in rownames(data_DE)) {
  x=ifelse(as.numeric(data_DE[gene,])!=0,yes = 1,no = 0)
  model_logistic=glm(x~lib_size_DE+condition_DE+batch_DE,family = "binomial")
  u=summary(model_logistic)
  u=c(u$coefficients[2:3,4],u$coefficients[2:3,1])
  model_selection=rbind(model_selection,u)
}

rownames(model_selection)=rownames(data_DE)
model_selection=data.frame(model_selection)
colnames(model_selection)=c("Lib_size","Condition","Intercep","Effect")
model_selection$Condition=p.adjust(model_selection$Condition,method = "BH")
model_selection$Lib_size=p.adjust(model_selection$Lib_size,method = "BH")
model_selection$Condition=-log10(model_selection$Condition)
model_selection$Lib_size=-log10(model_selection$Lib_size)
model_selection=model_selection[order(model_selection$Condition,decreasing = T),]

mean_DC_1_expression=rowMeans(TPM_data[,r$clusters$PCA$community==8])
mean_DC_2_expression=rowMeans(TPM_data[,r$clusters$PCA$community==9])
fit_DC_logFC=loess.smooth(y = mean_DC_2_expression,x = mean_DC_1_expression,evaluation = 1000,degree = 2,span = 1/20)
corrected_DC=c()
for (k in names(mean_DC_1_expression)) {
  u=which.min(abs(fit_DC_logFC$x-(mean_DC_1_expression[k])))
  corrected_DC=c(corrected_DC,mean_DC_2_expression[k]-fit_DC_logFC$y[u])
}

DC_valid_genes=rownames(model_selection)[model_selection$Lib_size>1.3]

DC_changed_genes=rownames(model_selection)[(model_selection$Condition>2 & model_selection$Lib_size>2)]
DC_changed_genes=intersect(DC_changed_genes,names(which(abs(corrected_DC)>2)))

plot(corrected_DC[DC_valid_genes],model_selection[DC_valid_genes,"Condition"],pch=16,
     xlab="DCs Mig. (2) vs DCs Mig. (1) (Log2 Fold Change)",ylab="-Log10(P-value)",cex.lab=1.5,
     col=string.to.colors(DC_valid_genes%in%DC_changed_genes,colors = c(alpha("red",0.7),alpha("black",0.3))))



data_expression=TPM_data[,r$clusters$PCA$community%in%order_cluster_bis]
cluster_expression=factor(r$clusters$PCA$community[r$clusters$PCA$community%in%order_cluster],order_cluster_bis)


####VI)Mean expression across clusters with violin plots

library(ggplot2)
violin_plot_gene=function(gene) {
  data_gene=data.frame(Expression=data_expression[gene,],
                       Condition=cluster_expression)
  ggplot(data_gene, aes(x=Condition, y=Expression))  +  geom_violin(trim=TRUE,scale = "width",fill='grey', color="black",bw=1.5)+ 
    scale_y_continuous(name = "Expression log2(TPM)",limits = c(0,max(data_expression[gene,]))) +
    theme_classic() + ggtitle(gene) + theme(plot.title = element_text(size=22),axis.text=element_text(size = 15),axis.title = element_text(size = 12)) + scale_x_discrete(labels=1:10,name=" ")  
}
View(r$diffgenes$PCA$community$`3`)

marker_list=c("Cd40","Cd86","Dll4","Il12b","Il18","Ccl17")
pdf("Project_Weizmann/Revision_plot/Violin_plot_Agpos.pdf",width = 4.5,height = 3.5)
for (k in c(marker_list)) {
  print(violin_plot_gene(k))
}
dev.off()

#VII)UMAP plot with various labels

#A)COnstruction of the plot

pdf("Project_Weizmann/Revision_plot/UMAP_ag_pos.pdf",width = 8,height = 8)
cell_of_interest=r$clusters$PCA$community%in%order_cluster
umap_plot=umap(r$reductions$PCA[cell_of_interest,],n_neighbors = 20,metric = "cosine",verbose = T,spread = 6)

#B)Plot themselves : Batch/Condition/Cluster

condition_count_bis=strsplit(condition_count,split = "-rep")
condition_count_bis=unlist(lapply(condition_count_bis,FUN = function(x){x[1]}))

plot(umap_plot,pch=16,xaxt="n",yaxt="n",bty="n",
     xlab="UMAP 1",ylab="UMAP 2",cex=0.7,main = "Louvain clustering",
     col=alpha(string.to.colors(r$clusters$PCA$community[cell_of_interest]),alpha = 0.3))

color_condition=c("cornflowerblue","salmon","darkblue","darkred","palegreen","darkgreen")

plot(umap_plot,pch=16,xaxt="n",yaxt="n",bty="n",
     xlab="UMAP 1",ylab="UMAP 2",cex=1,main="Condition ",
     col=alpha(string.to.colors(condition_count_bis[cell_of_interest],color_condition),alpha = 0.3))


plot(umap_plot,pch=16,xaxt="n",yaxt="n",bty="n",
     xlab="UMAP 1",ylab="UMAP 2",cex=0.7,main="Sequencing batch",
     col=alpha(string.to.colors(batch_count[cell_of_interest]),alpha = 0.3))

dev.off()

#C)Plot with gene labelling

plot(umap_plot,pch=16,xaxt="n",yaxt="n",bty="n",
     xlab="UMAP 1",ylab="UMAP 2",cex=1.5,
     col=alpha(color_convertion(r$counts[cell_of_interest,"Irg1"]),alpha = 0.1))

##VIII)Analysing Index Sorting / Flow cytometry data (INX data)

rescale_function=function(x) {
  x=(x-min(x,na.rm = T))
  x=x/max(x,na.rm = T)
  x=x*5
  return(x)
}

INX_data=read.table("Project_Weizmann/Final_Antigen_pos/FACS_data_final.txt") ### Read the file containing the INX data for Ag pos cells

INX_data_heatmap=INX_data[colnames(heatmap_data),]
rownames(INX_data_heatmap)=colnames(heatmap_data)


INX_data_heatmap$FSC.A=rescale_function(INX_data_heatmap$FSC.A)
INX_data_heatmap$SSC.A=rescale_function(INX_data_heatmap$SSC.A)

##We create a generic function to plot the INX data as shown in heatmap

plot_FACS_marker=function(marker="FSC.A",col="red") {
  plot(INX_data_heatmap[,marker],pch=21,bg=col,xaxt="n",xlab="",)
  abline(v=which(diff(as.numeric(r$clusters$PCA$community[colnames(heatmap_data)]))!=0),lwd=2,lty=2)
}

plot_FACS_marker("FSC.A")
plot_FACS_marker("SSC.A")
plot_FACS_marker("Antigen")
plot_FACS_marker("MHCII")
plot_FACS_marker("CD11b")
plot_FACS_marker("CD11c")



