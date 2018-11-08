
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

data_final=read.table("Project_Weizmann/Final_sequencing/data_final.txt",row.names = 1,sep="\t")

annotation=read.table("Project_Weizmann/Final_sequencing/annotation.txt",row.names = 1,sep="\t")
annotation=as.character(annotation$x)
names(annotation)=colnames(data_final)

condition=read.table("Project_Weizmann/Final_sequencing/condition.txt",row.names = 1,sep="\t")
condition=as.character(condition$x)
names(condition)=colnames(data_final)

batch=read.table("Project_Weizmann/Final_sequencing/batch.txt",row.names = 1,sep="\t")
batch=as.character(batch$x)
names(batch)=colnames(data_final)

replicates=read.table("Project_Weizmann/Final_sequencing/raw_condition.txt",row.names = 1,sep="\t")
replicates=as.character(replicates$x)
names(replicates)=colnames(data_final)

#####II)Filtering of the data

#A)Cell filtering -> removing cells with low amount of UMIs
lib_size=colSums(data_final)
names(lib_size)=colnames(data_final)
par(las=1)
hist(log10(1+lib_size),xlab="Library size of the cells (log10)",main="Librarys size distribution",n=100,xlim=c(1,8))
abline(v=log10(350),lwd=2,col="red",lty=2)

##B)Gene filtering  -> removing genes with low amount of UMIs

gene_size=rowSums(data_final)
hist(log10(1+gene_size),xlab="Gene abundance (log10)",main="Gene abundance distribution",n=100)
abline(v=log10(200),lwd=2,col="red",lty=2)
dev.off()

##C)Creation of the final dataset 

data_count=data_final[gene_size>250,lib_size>350] ###Pagoda2 only accepts Sparse Matrix as input
data_count=as(as.matrix(data_count),"dgCMatrix")
condition_count=condition[colnames(data_count)]
replicates_count=replicates[colnames(data_count)]
batch_count=batch[colnames(data_count)]

#III)Analysis using Pagoda2 pipeline

#A)Analysis per se

r <- Pagoda2$new(data_count,log.scale=FALSE)
r$adjustVariance(plot=T,gam.k=10)
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3)
r$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine')

r$getKnnClusters(method=multilevel.community,type='PCA')
r$getKnnClusters(method=infomap.community,type='PCA',name="infomap")
r$getKnnClusters(method=walktrap.community,type='PCA',name="walktrap")

#B)Visualisation using UMAP : UMAP is not implemented by PAGAODA2 pipeline intitially, R implementation by jlmelville

UMAP_all_cells=umap(r$reductions$PCA,spread = 6,
                        n_neighbors = 40,verbose =T,metric = "cosine")
##Spread parameter increased to improve visualisation


##Looking at clustering results projected on UMAP embedding 

plot(UMAP_all_cells,pch=16,xaxt="n",yaxt="n",bty="n",cex=0.7,
     xlab="UMAP 1",ylab="UMAP 2",main="Louvain clustering",
     col=alpha(string.to.colors(r$clusters$PCA$community),alpha = 0.4))

#IV)Creation of the heatmap (Figure 1B)

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
barplot(odd_score[1:10],horiz = T,xlim=c(0,0.5),main="Top 10 outlier score") ###

selected_odd_genes=names(which(odd_score>(mean(odd_score)+0.7*sd(odd_score))))

#B)Selecting marker genes 

r$getDifferentialGenes(type='PCA',verbose=T,clusterType='community',upregulated.only = T,z.threshold = 2)

gene_list=c()
gene_data_frame=c()

size_clusters=table(r$clusters$PCA$community)
order_cluster=order(size_clusters,decreasing=T)
size_clusters=(size_clusters)[order_cluster]
size_clusters=size_clusters/ncol(data_count)
order_cluster=order_cluster[size_clusters>0.01] ## We only taker clusters whose size represent at least more than 1% of total cell population
order_cluster=order_cluster[!order_cluster==16 & !order_cluster==7] 
## Removing the cells with the highest amount of Spike in : Low quality cells
## Also removing dividing cells : most of the UMIs are taken by dividing genes


for (k in order_cluster) {
  l =  (r$diffgenes$PCA$community)[[k]]
  l=l[l$highest==T,] ##ERCC and non-annotated genes are removed from marker list 
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
  
  l=l[intersect(rownames(l),selected_odd_genes),] ## Only genes with an excess of zeros are selected
  l=l[order(l$M,decreasing = T),] ###Genes are ranked based on their corresponding logFC compared to other cells
  if (nrow(l)==0) {
    next
  }
  l=l[selected_odd_genes,]
  gene_list=c(gene_list,rownames(l[1:5,])) ###
  gene_data_frame=rbind(gene_data_frame,data.frame(Gene=rownames(l),Score=l$Z,Cluster=k))
}
gene_list=unique(gene_list)
gene_list=intersect(gene_list,rownames(data_count))
gene_data_frame=na.omit(gene_data_frame)


#C)Computing normalised expression and concatenating gene expression into heatmap

TPM_data=t(log2(1+t(as.matrix(data_count))/lib_size[colnames(data_count)]*10^6))
TPM_data=as.matrix(TPM_data)

#Ordering clusters based on the size of the clusters

heatmap_data=c()
for (k in order_cluster) {
  v=TPM_data[gene_list,r$clusters$PCA$community==k]
  v=v[,sample(1:ncol(v),replace = F,size = ncol(v))]
  heatmap_data=cbind(heatmap_data,v)
  cat(paste(k,"\n"))
}

#D)Plot the heatmap

heatmap_data=as.matrix(heatmap_data)
lim_values=quantile(as.numeric((heatmap_data)),c(0,0.99)) ##Trimming extreme values

heatmap_data[heatmap_data<lim_values[1]]=lim_values[1]
heatmap_data[heatmap_data>lim_values[2]]=lim_values[2]

png("Project_Weizmann/Revision_plot/Heatmap_non_spec.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(130,40,4,4))
image(rotate(heatmap_data[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(order_cluster)])/sum(size_clusters[as.character(order_cluster)]),lwd=10,lty=2)
axis(2, at=seq(0,1,length.out=nrow((heatmap_data[,]))), 
     labels= ( rownames(heatmap_data)[nrow(heatmap_data):1]),
     las= 2, cex.axis=7, tick = F )
dev.off()


#V)Analysis of the proportion of cells across cell

#A) Computing major cell type proportion across the samples of intereste (control samples + day1 + day2)

##We want to select only the most important populations

size_clusters=table(r$clusters$PCA$community)
order_cluster=order(size_clusters,decreasing=T)
size_clusters=(size_clusters)[order_cluster]
size_clusters=size_clusters/ncol(data_count)
order_cluster=order_cluster[size_clusters>0.01]
order_cluster=order_cluster[!order_cluster==16 & !order_cluster==7]


cell_proportion=table(replicates_count,r$clusters$PCA$community)
cell_proportion=cell_proportion[c(7,6,9,10,11,13,16:17,19,20,1:4),order_cluster[-length(order_cluster)]]
cell_proportion=cbind(cell_proportion[,1:7],rowSums(cell_proportion[,8:9])) ### We merge all the small cluster : Macrophages, Plasma cells, Neutrophils
cell_proportion=cell_proportion/rowSums(cell_proportion)*100

rownames(cell_proportion)=c("PBS rep.1","PBS rep.2",
                            "Ms day1 rep.1","Ms day1 rep.2",
                            "Ms day2 rep.1","Ms day2 rep.2",
                            "Nb day1 rep.1","Nb day1 rep.2",
                            "Nb day2 rep.1","Nb day2 rep.2",
                            "Ca day1 rep.1","Ca day1 rep.2",
                            "Ca day2 rep.1","Ca day2 rep.2")

#B)Plot the proportion using stacked barplots

tick_barplot=barplot(t(cell_proportion[1:6,]),space = c(0.3,0.3,1,0.3,1,0.3),
                     width = 12,xlim=c(0,100),ylab="Cell type abundance (%)",xaxt="n",density = (8:1)*10)
text(cex=1, y = -10,x = tick_barplot-10, rownames(cell_proportion)[1:6], xpd=TRUE, srt=45)

tick_barplot=barplot(t(cell_proportion[c(1:2,7:10),]),space = c(0.3,0.3,1,0.3,1,0.3),
                     width = 12,xlim=c(0,100),ylab="Cell type abundance (%)",xaxt="n",density = (8:1)*10)
text(cex=1, y = -10,x = tick_barplot-10, rownames(cell_proportion[c(1:2,7:10),]), xpd=TRUE, srt=45)


tick_barplot=barplot(t(cell_proportion[c(1:2,11:14),]),space = c(0.3,0.3,1,0.3,1,0.3),
                     width = 12,xlim=c(0,100),ylab="Cell type abundance (%)",xaxt="n",density = (8:1)*10)
text(cex=1, y = -10,x = tick_barplot-10, rownames(cell_proportion[c(1:2,11:14),]), xpd=TRUE, srt=45)


#VI)Reproductbility of the study

#A)Construction of the datasets : TPM expression in NK, monocytes and Migratory DCs

monocytes_expression=aggregate(t(TPM_data[,r$clusters$PCA$community==13]),by=list(replicates_count[r$clusters$PCA$community==13]),FUN=mean)
monocytes_expression=monocytes_expression[c(7,8,9,10,11,13,16:17,19,20,1:4),]
rownames(monocytes_expression)=monocytes_expression[,1]
monocytes_expression=monocytes_expression[,-1]
monocytes_expression=t(monocytes_expression)

NK_expression=aggregate(t(TPM_data[,r$clusters$PCA$community==3]),by=list(replicates_count[r$clusters$PCA$community==3]),FUN=mean)
NK_expression=NK_expression[c(7,8,9,10,11,13,16:17,19,20,1:4),]
rownames(NK_expression)=NK_expression[,1]
NK_expression=NK_expression[,-1]
NK_expression=t(NK_expression)

DC_expression=aggregate(t(TPM_data[,r$clusters$PCA$community==6]),by=list(replicates_count[r$clusters$PCA$community==6]),FUN=mean)
DC_expression=DC_expression[c(7,8,9,10,11,13,16:17,19,20,1:4),]
rownames(DC_expression)=DC_expression[,1]
DC_expression=DC_expression[,-1]
DC_expression=t(DC_expression)

#B)Plot of the resulsts

par(las=1)
plot(monocytes_expression[,5:6],pch=16,cex=1,
     xlab="Ms day2 rep.1",ylab="Ms day2 rep.2",
     xlim=c(0,16),ylim=c(0,16),main = "Monocytes",col=alpha("black",0.2))
R=cor(monocytes_expression[,5:6])[2,1]
legend("topleft",legend = paste("R = ",round(R,digits = 2),sep = ""),bty="n")
abline(lm(monocytes_expression[,6]~monocytes_expression[,5]),lwd=2,lty=2)

plot(NK_expression[,5:6],pch=16,cex=1,
     xlab="Ms day2 rep.1",ylab="Ms day2 rep.2",
     xlim=c(0,16),ylim=c(0,16),main = "NK",col=alpha("black",0.2))
R=cor(NK_expression[,5:6])[2,1]
legend("topleft",legend = paste("R = ",round(R,digits = 2),sep = ""),bty="n")
abline(lm(NK_expression[,6]~NK_expression[,5]),lwd=2,lty=2)

plot(DC_expression[,5:6],pch=16,cex=1,
     xlab="Ms day2 rep.1",ylab="Ms day2 rep.2",
     xlim=c(0,16),ylim=c(0,16),main = "DCs Mig",col=alpha("black",0.2))
abline(lm(DC_expression[,6]~DC_expression[,5]),lwd=2,lty=2)
R=cor(DC_expression[,5:6])[2,1]
legend("topleft",legend = paste("R = ",round(R,digits = 2),sep = ""),bty="n")
dev.off()

#VII)Effects of the pathogens on the cells

#A)NK cells : DE genes detection

#1)Corrected logFCs

mean_expression_NK=aggregate(t(TPM_data[,r$clusters$PCA$community==6]),
                             by=list(factor(condition_count[r$clusters$PCA$community==6])),FUN=mean)
rownames(mean_expression_NK)=mean_expression_NK[,1]
mean_expression_NK=t(mean_expression_NK[,-1])
mean_expression_NK=mean_expression_NK[,c(4,6)]
plot(mean_expression_NK)
NK_loess=loess.smooth(y = mean_expression_NK[,2],x = mean_expression_NK[,1],degree = 2,span = 0.3,evaluation = 1000)
lines(NK_loess,lwd=2,col="red",lty=2)
abline(0,1)

corrected_NK=c()
for (k in rownames(mean_expression_NK)) {
  u=which.min(abs(NK_loess$x-(mean_expression_NK[k,1])))
  corrected_NK=c(corrected_NK,mean_expression_NK[k,2]-NK_loess$y[u])
}
names(corrected_NK)=rownames(mean_expression_NK)

#2)Logistic test

data_DE=data_count[,r$clusters$PCA$community==6 & condition_count %in% c("Control","M.smeg day2")]
condition_DE=condition_count[colnames(data_DE)]
lib_size_DE=lib_size[colnames(data_DE)]

model_selection_NK=c()
for (gene in rownames(data_DE)) {
  x=ifelse(as.numeric(data_DE[gene,])!=0,yes = 1,no = 0)
  model_logistic=glm(x~lib_size_DE+condition_DE,family = "binomial")
  u=summary(model_logistic)
  u=c(u$coefficients[2:3,4],u$coefficients[2:3,1])
  model_selection_NK=rbind(model_selection_NK,u)
}

rownames(model_selection_NK)=rownames(data_DE)
model_selection_NK=data.frame(model_selection_NK)
colnames(model_selection_NK)=c("Lib_size","Condition","Intercep","Effect")
model_selection_NK$Condition=p.adjust(model_selection_NK$Condition,method = "BH")
model_selection_NK$Lib_size=p.adjust(model_selection_NK$Lib_size,method = "BH")
model_selection_NK$Condition=-log10(model_selection_NK$Condition)
model_selection_NK$Lib_size=-log10(model_selection_NK$Lib_size)
model_selection_NK=model_selection_NK[order(model_selection_NK$Condition,decreasing = T),]

model_selection_NK_bis=model_selection_NK
model_selection_NK_bis=model_selection_NK_bis[model_selection_NK_bis$Lib_size>2 & !grepl("ERCC",rownames(model_selection_NK_bis)),]

u=model_selection_NK
par(las=1)
plot(corrected_NK[rownames(model_selection_NK_bis)],(model_selection_NK_bis$Condition),
     pch=16,xlim=c(-3.5,3.5),xlab="Ms day2 vs PBS (Log2 Fold Change)",ylab="-Log10(P-value)",cex=1,cex.lab=1.5,
     col=string.to.colors(model_selection_NK_bis$Condition>2,colors = c(alpha("red",alpha = 0.7),alpha("black",alpha = 0.3))))
grid(nx = 4,ny=4,lty=2)
plot(corrected_NK[rownames(model_selection_NK_bis)],(model_selection_NK_bis$Condition),
     pch=16,xlim=c(-3.5,3.5),xlab="Ms day2 vs PBS (Log2 Fold Change)",ylab="-Log10(P-value)",cex=1,
     col=alpha(string.to.colors(model_selection_NK_bis$Condition>2,colors = c("red","black")),alpha = 0.3))
text(corrected_NK[rownames(model_selection_NK_bis)],(model_selection_NK_bis$Condition),labels = rownames(model_selection_NK_bis) )

#B)Monocytes : looking for DE genes

#1)Corrected logFCs

mean_expression_mono=aggregate(t(TPM_data[,r$clusters$PCA$community==13]),
                               by=list(factor(condition_count[r$clusters$PCA$community==13])),FUN=mean)
rownames(mean_expression_mono)=mean_expression_mono[,1]
mean_expression_mono=t(mean_expression_mono[,-1])
mean_expression_mono=mean_expression_mono[,c(4,6)]
plot(mean_expression_mono)
mono_loess=loess.smooth(y = mean_expression_mono[,2],x = mean_expression_mono[,1],degree = 2,span = 0.1,evaluation = 1000)
lines(mono_loess,lwd=2,col="red",lty=2)
abline(0,1)

corrected_mono=c()
for (k in rownames(mean_expression_mono)) {
  u=which.min(abs(mono_loess$x-(mean_expression_mono[k,1])))
  corrected_mono=c(corrected_mono,mean_expression_mono[k,2]-mono_loess$y[u])
}
names(corrected_mono)=rownames(mean_expression_mono)

#2)Logistic test

data_DE=data_count[,r$clusters$PCA$community==13 & condition_count %in% c("Control","Nb day2")]
condition_DE=condition_count[colnames(data_DE)]
lib_size_DE=lib_size[colnames(data_DE)]

model_selection_mono=c()
for (gene in rownames(data_DE)) {
  x=ifelse(as.numeric(data_DE[gene,])!=0,yes = 1,no = 0)
  model_logistic=glm(x~lib_size_DE+condition_DE,family = "binomial")
  u=summary(model_logistic)
  u=c(u$coefficients[2:3,4],u$coefficients[2:3,1])
  model_selection_mono=rbind(model_selection_mono,u)
}

rownames(model_selection_mono)=rownames(data_DE)
model_selection_mono=data.frame(model_selection_mono)
colnames(model_selection_mono)=c("Lib_size","Condition","Intercep","Effect")
model_selection_mono$Condition=p.adjust(model_selection_mono$Condition,method = "BH")
model_selection_mono$Lib_size=p.adjust(model_selection_mono$Lib_size,method = "BH")
model_selection_mono$Condition=-log10(model_selection_mono$Condition)
model_selection_mono$Lib_size=-log10(model_selection_mono$Lib_size)
model_selection_mono=model_selection_mono[order(model_selection_mono$Condition,decreasing = T),]

model_selection_mono_bis=model_selection_mono
model_selection_mono_bis=model_selection_mono_bis[model_selection_mono_bis$Lib_size>2 & !grepl("ERCC",rownames(model_selection_mono_bis)),]

par(las=1)
plot(corrected_mono[rownames(model_selection_mono_bis)],(model_selection_mono_bis$Condition),
     pch=16,xlab="Nb day2 vs PBS (Log2 Fold Change)",ylab="-Log10(P-value)",cex=1,xlim=c(-9,9),cex.lab=1.5,
     col=string.to.colors(model_selection_mono_bis$Condition>2,colors = c(alpha("red",alpha = 0.7),alpha("black",alpha = 0.3))))
grid(nx = 4,ny=4,lty=2)
plot(corrected_mono[rownames(model_selection_mono_bis)],(model_selection_mono_bis$Condition),
     pch=16,xlab="Log2FC (Ms day2 vs PBS)",ylab="-Log10(P-value)",cex=1,
     col=string.to.colors(model_selection_mono_bis$Condition>2,colors = c(alpha("red",alpha = 0.7),alpha("black",alpha = 0.3))))
text(corrected_mono[rownames(model_selection_mono_bis)],(model_selection_mono_bis$Condition),labels =rownames(model_selection_mono_bis))


#VIII)Effects of the pathogens on the cellular composition

#A)Global description

cell_type_count=table(replicates_count,r$clusters$PCA$community)
cell_type_count=cell_type_count[c(7,6,9,10,11,13,16:17,19,20,1:4),order_cluster]
cell_type_count=cbind(cell_type_count[,1:7],rowSums(cell_type_count[,8:10])) ### We merge all the small clusters : Macrophages, Plasma cells, Neutrophils
colnames(cell_type_count)=c("Monocytes","MigDCs","NK","ResDCs 2","Lymphocytes",
                            "pDCs","ResDCs 1","Other cells")
total_cell=rowSums(cell_type_count)

estimated_p_parameters=cell_type_count[1:2,]
estimated_p_parameters=estimated_p_parameters/rowSums(estimated_p_parameters)
estimated_p_parameters=colMeans(estimated_p_parameters)

#B)Statistical significance 


binom_matrix=matrix(0,nrow = nrow(cell_abundance),ncol=ncol(cell_abundance) )

for (i in 1:nrow(cell_type_count)) {
  for (j in 1:ncol(cell_type_count)) {
    u=binom.test(x = cell_type_count[i,j],p =estimated_p_parameters[j],n = total_cell[j],alternative="greater")
    
    binom_matrix[i,j]=u$p.value
  }
}

binom_matrix=matrix(p.adjust(binom_matrix,method = "bonf"),nrow = nrow(binom_matrix),ncol = ncol(binom_matrix) )
binom_matrix=-log10(binom_matrix)
binom_matrix=round(binom_matrix,digits = 2)

#C)Enrichment analysis 

cell_type_proportion=round(cell_type_count/total_cell,2)
cell_type_enrichment=t(round(t(cell_type_proportion)/estimated_p_parameters,digits = 2))

rownames(cell_type_proportion)=c("PBS rep.1","PBS rep.2",
                                 "Ms day1 rep.1","Ms day1 rep.2",
                                 "Ms day2 rep.1","Ms day2 rep.2",
                                 "Nb day1 rep.1","Nb day1 rep.2",
                                 "Nb day2 rep.1","Nb day2 rep.2",
                                 "Ca day1 rep.1","Ca day1 rep.2",
                                 "Ca day2 rep.1","Ca day2 rep.2")
colnames(cell_type_proportion)=c("Monocytes","DCs Mig.","NK","DC2s Res.","Lymphocytes","pDCs","DC1s Res.","Others")

colnames(binom_matrix)=colnames(cell_type_proportion)
rownames(binom_matrix)=rownames(cell_type_proportion)

colnames(cell_type_enrichment)=colnames(cell_type_proportion)
rownames(cell_type_enrichment)=rownames(cell_type_proportion)

#D)Bubble plot
 
x_position=1:ncol(cell_type_proportion)+2
y_position=nrow(cell_type_proportion):1
pairwise_position=expand.grid(x_position,y_position)
par(mar=c(2,2,2,2),family="Arial")
plot(pairwise_position,cex=as.numeric((as.matrix(t(cell_type_proportion))))*8+0.5, pch=16,
     col=color_convertion(t(binom_matrix),max_scale = quantile(binom_matrix,probs = 0.9)),
     bty="n",xaxt="n",yaxt="n",xlim=c(0,ncol(cell_type_proportion)+1),ylim=c(-1.5,nrow(cell_type_proportion)+2.5),xlab="",ylab="")
points(pairwise_position,cex=as.numeric((as.matrix(t(cell_type_proportion))))*8+0.5,lwd=1.4)
abline(h=c(2.5,4.5,6.5,8.5,10.5,12.5),lwd=2,lty=2)
points(x = c(1,2,3.2),y=rep(nrow(cell_type_proportion)+2,3),cex=c(0.05,0.2,0.4)*8+0.5,lwd=1.4)
text(x = c(1,2,3.2),y=rep(nrow(cell_type_proportion)+3,3),labels = c("5%","20%","40%"))
text(x = 1,y = y_position,labels = rownames(cell_type_proportion))
text(x = x_position-0.9,y = -0.8,labels = colnames(cell_type_proportion),srt=45,cex=1.2)
color.legend(xl = 5,yb = nrow(cell_type_proportion)+2.3,xr = 9,yt = nrow(cell_type_proportion)+1.6,legend="",
             rect.col = colorRampPalette(c("white","darkorange"))(100))
text(x = ncol(cell_type_proportion)-1,y=nrow(cell_type_proportion)+2.7, labels = "-Log10(p-value)")


####VII)Mean expression across clusters with violin plots

library(ggplot2)

data_expression=TPM_data[,r$clusters$PCA$community%in%order_cluster]
cluster_expression=factor(r$clusters$PCA$community[r$clusters$PCA$community%in%order_cluster],order_cluster)

violin_plot_gene=function(gene) {
  data_gene=data.frame(Expression=data_expression[gene,],
                       Condition=cluster_expression)
  ggplot(data_gene, aes(x=Condition, y=Expression))  +  geom_violin(trim=TRUE,scale = "width",fill='grey', color="black",bw=1.5)+ 
    scale_y_continuous(name = "Expression log2(TPM)",limits = c(0,max(data_expression[gene,]))) +
    theme_classic() + ggtitle(gene) + theme(plot.title = element_text(size=22),axis.text=element_text(size = 15),axis.title = element_text(size = 12)) + scale_x_discrete(labels=1:10,name=" ")  
}
par(las=1,mfrow=c(4,2))

extra_marker_list_2=c("Fcgr2b","Ccl22","Gzma","Trbc2","Siglech","Xcr1","S100a9","Apoe")

pdf("Project_Weizmann/Revision_plot/Violin_plot_non_spec.pdf",width = 4.5,height = 3.5)
for (k in c(extra_marker_list_2)) {
  print(violin_plot_gene(k))
}
dev.off()

#IX)UMAP plot with various labels

#A)Creation of the projection : only project cells from the 10 clusters and 7 conditions studied

replicate_list=levels(factor(replicates_count))
replicate_list=replicate_list[c(7,6,9,10,11,13,16:17,19,20,1:4)]
cell_of_interest= (r$clusters$PCA$community %in% order_cluster) & (replicates_count %in%replicate_list )

UMAP_cell_interest=umap(r$reductions$PCA[cell_of_interest,],spread = 6,
                        n_neighbors = 40,verbose =T,metric = "cosine")

#B)Ploting the batch/conditions/clusters

plot(UMAP_cell_interest,pch=16,xaxt="n",yaxt="n",bty="n",cex=0.7,
     xlab="UMAP 1",ylab="UMAP 2",main="Louvain clustering",
     col=alpha(string.to.colors(r$clusters$PCA$community[cell_of_interest]),alpha = 0.4))

plot(UMAP_cell_interest,pch=16,xaxt="n",yaxt="n",bty="n",cex=0.7,
     xlab="UMAP 1",ylab="UMAP 2",main="Infomap clustering",
     col=alpha(string.to.colors(r$clusters$PCA$infomap[cell_of_interest]==3),alpha = 0.4))


plot(UMAP_cell_interest,pch=16,xaxt="n",yaxt="n",bty="n",cex=0.7,
     xlab="UMAP 1",ylab="UMAP 2",main="Sequencing batch",
     col=alpha(string.to.colors(batch_count[cell_of_interest]),alpha = 0.4))

color_condition=c("darkgreen","darkred","grey","cornflowerblue","darkblue","palegreen","salmon")

plot(UMAP_cell_interest,pch=16,xaxt="n",yaxt="n",bty="n",cex=0.7,
     xlab="UMAP 1",ylab="UMAP 2",main="Conditions",
     col=alpha(string.to.colors(condition_count[cell_of_interest],color_condition),alpha = 0.4))

legend("bottomleft",legend = c("PBS",
                               "Ms day1","Ms day2",
                               "Nb day1","Nb day2",
                               "Ca day1","Ca day2"),
       col = c("grey","salmon","darkred",
               "cornflowerblue","darkblue",
               "palegreen","darkgreen"),bty="n",pch=16)

#C)Ploting one specificic genes

plot(UMAP_cell_interest,pch=16,xaxt="n",yaxt="n",bty="n",cex=0.7,
     xlab="UMAP 1",ylab="UMAP 2",
     col=alpha(color_convertion(r$counts[cell_of_interest,"Bst2"]),alpha = 0.4))

##X)Analysing Index Sorting / Flow cytometry data (INX data)

rescale_function=function(x) {
  x=(x-min(x,na.rm = T))
  x=x/max(x,na.rm = T)
  x=x*5
  return(x)
}

INX_data=read.table("Project_Weizmann/Final_sequencing/Pagoda_2_results/INX_data.txt") ### Read the file containing the INX data for total cell data

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



