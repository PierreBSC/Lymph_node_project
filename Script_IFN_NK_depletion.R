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

data_final=read.table("Project_Weizmann/NK_depletion_validation/data_final.txt",row.names = 1,sep="\t")

condition=read.table("Project_Weizmann/NK_depletion_validation/condition.txt",row.names = 1,sep="\t")
condition=as.character(condition$x)
names(condition)=colnames(data_final)

replicates=read.table("Project_Weizmann/NK_depletion_validation/replicates.txt",row.names = 1,sep="\t")
replicates=as.character(replicates$x)
names(replicates)=colnames(data_final)

annotation=read.table("Project_Weizmann/NK_depletion_validation/annotation.txt",row.names = 1,sep="\t")
annotation=as.character(annotation$x)
names(annotation)=colnames(data_final)

##### II) Analysis of the Ag+ fraction

data_ag_pos=data_final[,annotation=="Antigen+"]
condition_ag_pos=condition[annotation=="Antigen+"]
replicates_ag_pos=replicates[annotation=="Antigen+"]

#A)QC and filtering

lib_size=colSums(data_ag_pos)
names(lib_size)=colnames(data_ag_pos)
hist()

par(las=1,family="serif")
hist(log10(1+lib_size),xlab="Library size of the cells (log10)",main="Librarys size distribution",n=100,xlim=c(1,8))
abline(v=log10(500),lwd=1,col="black",lty=2)
abline(v=log10(350),lwd=2,col="red",lty=2)


gene_size=rowSums(data_ag_pos)
hist(log10(1+gene_size),xlab="Gene abundance (log10)",main="Gene abundance distribution",n=100)
abline(v=log10(50),lwd=2,col="red",lty=2)


data_ag_pos_count=data_ag_pos[gene_size>50,lib_size>350]
data_ag_pos_count=as(as.matrix(data_ag_pos_count),"dgCMatrix")
condition_ag_pos_count=condition_ag_pos[colnames(data_ag_pos_count)]
replicates_ag_pos_count=replicates_ag_pos[colnames(data_ag_pos_count)]

#B)Analysis per se 

r <- Pagoda2$new(data_ag_pos_count,log.scale=FALSE)
r$adjustVariance(plot=T,gam.k=10)
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3)
r$makeKnnGraph(k=20,type='PCA',center=T,distance='cosine')

r$getKnnClusters(method=multilevel.community,type='PCA')
r$getKnnClusters(method=infomap.community,type='PCA',name="infomap",)
r$getKnnClusters(method=walktrap.community,type='PCA',name="walktrap")

#C)Heatmap creation

proportion_zero=rowSums(as.matrix(data_ag_pos_count)==0)/ncol(data_ag_pos_count)
mean_expression=rowMeans(as.matrix(data_ag_pos_count))

par(las=1,family="serif")
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
barplot(odd_score[1:10],horiz = T,xlim=c(0,0.5),main="Top 10 outlier score")

selected_odd_genes=names(which(odd_score>(mean(odd_score)+2*sd(odd_score))))

r$getDifferentialGenes(type='PCA',verbose=T,clusterType='infomap',upregulated.only = T,z.threshold = 3)

gene_list=c()
gene_data_frame=c()

size_clusters=table(r$clusters$PCA$infomap)
order_cluster=order(size_clusters,decreasing=T)
size_clusters=(size_clusters)[order_cluster]
size_clusters=size_clusters/ncol(data_non_spec_count)
order_cluster=order_cluster[size_clusters>0.01]
order_cluster=order_cluster[order_cluster!=2 & order_cluster!=6 ] ##Removing higly dividing cells and ERCC high cluster


for (k in order_cluster) {
  l =  (r$diffgenes$PCA$info)[[k]]
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
  
  l=l[intersect(rownames(l),selected_odd_genes),]
  if (nrow(l)==0) {
    next
  }
  l=l[order(l$M,decreasing = T),]
  
  gene_list=c(gene_list,rownames(l[1:5,]))
  gene_data_frame=rbind(gene_data_frame,data.frame(Gene=rownames(l),Score=l$Z,Cluster=k))
}
gene_list=unique(gene_list)
gene_list=intersect(gene_list,rownames(data_ag_pos_count))
gene_data_frame=na.omit(gene_data_frame)

#r$plotGeneHeatmap(genes = gene_list,type = "PCA",clusterType = "infomap",cluster.genes = F)
lib_size=colSums(data_ag_pos)
TPM_data=t(log2(1+t(as.matrix(data_ag_pos_count))/lib_size[colnames(data_ag_pos_count)]*10^6))
TPM_data=as.matrix(TPM_data)

heatmap_data_ag_pos=c()
for (k in order_cluster) {
  v=TPM_data[intersect(gene_list,selected_odd_genes),r$clusters$PCA$infomap==k]
  v=v[,sample(1:ncol(v),replace = F,size = ncol(v))]
  heatmap_data_ag_pos=cbind(heatmap_data_ag_pos,v)
  cat(paste(k,"\n"))
}

###Plot the heatmap

heatmap_data_ag_pos=as.matrix(heatmap_data_ag_pos)
lim_values=quantile(as.numeric((heatmap_data_ag_pos)),c(0.0,0.99))

heatmap_data_ag_pos[heatmap_data_ag_pos<lim_values[1]]=lim_values[1]
heatmap_data_ag_pos[heatmap_data_ag_pos>lim_values[2]]=lim_values[2]

png("Project_Weizmann/Revision_plot/Heatmap_NK_dep_ag_pos.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(4,40,320,4))
image(rotate(heatmap_data_ag_pos[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(order_cluster)])/sum(size_clusters[as.character(order_cluster)]),lwd=10,lty=2)
axis(2, at=seq(0,1,length.out=nrow((heatmap_data_ag_pos[,]))), 
     labels= ( rownames(heatmap_data_ag_pos)[nrow(heatmap_data_ag_pos):1]),
     las= 2, cex.axis=5, tick = F )
dev.off()


ag_pos_proportion=table(replicates_ag_pos_count[colnames(heatmap_data_ag_pos)],
                        as.numeric(as.character(r$clusters$PCA$infomap[colnames(heatmap_data_ag_pos)])))
ag_pos_proportion=ag_pos_proportion/rowSums(ag_pos_proportion)
black_palette=colorRampPalette(colors = c("black","white"))(4)
barplot(t(ag_pos_proportion),col=black_palette)

##### III) Analysis of the non B/T fraction

data_non_spec=data_final[,annotation=="non B/T"]
condition_non_spec=condition[annotation=="non B/T"]
replicates_non_spec=replicates[annotation=="non B/T"]

#A)QC and filtering

lib_size=colSums(data_non_spec)
names(lib_size)=colnames(data_non_spec)
par(las=1,family="serif")
hist(log10(1+lib_size),xlab="Library size of the cells (log10)",main="Librarys size distribution",n=100,xlim=c(1,5))
abline(v=log10(500),lwd=1,col="black",lty=2)
abline(v=log10(350),lwd=2,col="red",lty=2)


gene_size=rowSums(data_non_spec)
hist(log10(1+gene_size),xlab="Gene abundance (log10)",main="Gene abundance distribution",n=100)
abline(v=log10(50),lwd=2,col="red",lty=2)

data_non_spec_count=data_non_spec[gene_size>50,lib_size>350]
data_non_spec_count=as(as.matrix(data_non_spec_count),"dgCMatrix")
condition_non_spec_count=condition_non_spec[colnames(data_non_spec_count)]
replicates_non_spec_count=replicates_non_spec[colnames(data_non_spec_count)]

#B)Analysis per se 

r_bis <- Pagoda2$new(data_non_spec_count,log.scale=FALSE)
r_bis$adjustVariance(plot=T,gam.k=10)
r_bis$calculatePcaReduction(nPcs=100,n.odgenes=3e3)
r_bis$makeKnnGraph(k=20,type='PCA',center=T,distance='cosine')

r_bis$getKnnClusters(method=multilevel.community,type='PCA')
r_bis$getKnnClusters(method=infomap.community,type='PCA',name="infomap",)
r_bis$getKnnClusters(method=walktrap.community,type='PCA',name="walktrap")

#C)Heatmap creation

proportion_zero=rowSums(as.matrix(data_non_spec_count)==0)/ncol(data_non_spec_count)
mean_expression=rowMeans(as.matrix(data_non_spec_count))

par(las=1,family="serif")
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
barplot(odd_score[1:10],horiz = T,xlim=c(0,0.5),main="Top 10 outlier score")

selected_odd_genes=names(which(odd_score>(mean(odd_score)+0.5*sd(odd_score))))

r_bis$getDifferentialGenes(type='PCA',verbose=T,clusterType='infomap',upregulated.only = T,z.threshold = 2)

size_clusters=table(r_bis$clusters$PCA$infomap)
order_cluster=order(size_clusters,decreasing=T)
size_clusters=(size_clusters)[order_cluster]
size_clusters=size_clusters/ncol(data_non_spec_count)
order_cluster=order_cluster[size_clusters>0.05]
order_cluster=order_cluster[order_cluster!=5 & order_cluster!=10]
gene_list=c()
gene_data_frame=c()

for (k in order_cluster) {
  l =  (r_bis$diffgenes$PCA$infomap)[[k]]
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
gene_list=intersect(gene_list,rownames(data_ag_pos_count))
gene_data_frame=na.omit(gene_data_frame)

r_bis$plotGeneHeatmap(genes = gene_list,type = "PCA",clusterType = "infomap",cluster.genes = F)

lib_size=colSums(data_non_spec)
TPM_data=t(log2(1+t(as.matrix(data_non_spec_count))/lib_size[colnames(data_non_spec_count)]*10^6))
TPM_data=as.matrix(TPM_data)

heatmap_data_non_spec=c()
for (k in order_cluster) {
  v=TPM_data[intersect(gene_list,selected_odd_genes),r_bis$clusters$PCA$infomap==k]
  v=v[,sample(1:ncol(v),replace = F,size = ncol(v))]
  heatmap_data_non_spec=cbind(heatmap_data_non_spec,v)
  cat(paste(k,"\n"))
}


###Plot the heatmap

heatmap_data_non_spec=as.matrix(heatmap_data_non_spec)
lim_values=quantile(as.numeric((heatmap_data_non_spec)),c(0.0,0.99))

heatmap_data_non_spec[heatmap_data_non_spec<lim_values[1]]=lim_values[1]
heatmap_data_non_spec[heatmap_data_non_spec>lim_values[2]]=lim_values[2]

png("Project_Weizmann/Revision_plot/Heatmap_NK_dep_non_spec.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(200,40,4,4))
image(rotate(heatmap_data_non_spec[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
box(which = "plot",lty = 1,lwd=10)
abline(v=cumsum(size_clusters[as.character(order_cluster)])/sum(size_clusters[as.character(order_cluster)]),lwd=10,lty=2)
axis(2, at=seq(0,1,length.out=nrow((heatmap_data_non_spec[,]))), 
     labels= ( rownames(heatmap_data_non_spec)[nrow(heatmap_data_non_spec):1]),
     las= 2, cex.axis=7, tick = F )
dev.off()


non_spec_proportion=table(replicates_non_spec_count[colnames(heatmap_data_non_spec)],
                          as.numeric(as.character(r_bis$clusters$PCA$infomap[colnames(heatmap_data_non_spec)])))
non_spec_proportion=cbind(non_spec_proportion[,1:7],rowSums(non_spec_proportion[,8:10]))

non_spec_proportion=non_spec_proportion/rowSums(non_spec_proportion)
black_palette=colorRampPalette(colors = c("black","white"))(8)
barplot(t(non_spec_proportion),col=black_palette)


#####IV) Looking at the effect of the depletion on possible target cells


#A)Detction of DE genes using a logistic based appraoch 


#1)General function :

logistic_DE=function(cluster,gating="Ag_pos") {
  
  if (gating=="Ag_pos") {
    data_DE=data_ag_pos_count[,r$clusters$PCA$infomap==cluster]
    condition_DE=condition_ag_pos_count[r$clusters$PCA$infomap==cluster]
    condition_DE=factor(condition_DE,levels = c("Myco day2 + isotype control","Myco day2 + anti-Nk1.1/Ifng"))
    lib_size_DE=log10(colSums(data_final[,colnames(data_DE)]))
  }
  
  if (gating=="Non_spec") {
    data_DE=data_non_spec_count[,r_bis$clusters$PCA$infomap==cluster]
    condition_DE=condition_non_spec_count[r_bis$clusters$PCA$infomap==cluster]
    condition_DE=factor(condition_DE,levels = c("Myco day2 + isotype control","Myco day2 + anti-Nk1.1/Ifng"))
    lib_size_DE=log10(colSums(data_final[,colnames(data_DE)]))
  }
  
  
  model_selection=c()
  for (gene in rownames(data_DE)) {
    x=ifelse(as.numeric(data_DE[gene,])!=0,yes = 1,no = 0)
    model_logistic=glm(x~lib_size_DE+condition_DE,family = "binomial")
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
  cat(paste(round(sum(model_selection$Lib_size>2)/length(model_selection$Lib_size),digits = 2),"percents of the genes fit to the model"))
  return(list(table=model_selection,up_regulated_genes=rownames(model_selection[model_selection$Lib_size>3 & model_selection$Condition>1.3 & model_selection$Effect>0,]),
              down_regulated_genes=rownames(model_selection[model_selection$Lib_size>3 & model_selection$Condition>1.3 & model_selection$Effect<0,])))
}

#2)Looking at antigen positive cells

DE_macro_infla=logistic_DE(cluster = 1,gating = "Ag_pos")
DE_DC_ag_pos=logistic_DE(cluster = 3,gating = "Ag_pos")
DE_macro_resident=logistic_DE(cluster = 4,gating = "Ag_pos")
DE_neutro=logistic_DE(cluster = 5,gating = "Ag_pos")

#3)Looking at non-specifics cells

DE_mono=logistic_DE(cluster = 1,gating = "Non_spec")
DE_DC_non_spec=logistic_DE(cluster = 2,gating = "Non_spec")
DE_NK=logistic_DE(cluster = 3,gating = "Non_spec")
DE_ILC=logistic_DE(cluster = 4,gating = "Non_spec")
DE_DC_2=logistic_DE(cluster = 6,gating = "Non_spec")
DE_pDC=logistic_DE(cluster = 7,gating = "Non_spec")
DE_DC_3=logistic_DE(cluster = 8,gating = "Non_spec")

#B)Plots 

#1)Barplot to look at the total number of DE genes 

get_DE_genes_number=function(x) {
  up=sum((x$Lib_size>1.3 & x$Condition>1.3 & x$Effect>0))
  down=sum((x$Lib_size>1.3 & x$Condition>1.3 & x$Effect<0))
  return(c(down,up))
}

effect_ag_pos=cbind(get_DE_genes_number(DE_macro_infla$table),
                    get_DE_genes_number(DE_DC_ag_pos$table),
                    get_DE_genes_number(DE_macro_resident$table),
                    get_DE_genes_number(DE_neutro$table))

effect_non_spec=cbind(get_DE_genes_number(DE_mono$table),
                      get_DE_genes_number(DE_DC_non_spec$table),
                      get_DE_genes_number(DE_NK$table),
                      get_DE_genes_number(DE_ILC$table),
                      get_DE_genes_number(DE_DC_2$table),
                      get_DE_genes_number(DE_pDC$table),
                      get_DE_genes_number(DE_DC_3$table))

colnames(effect_ag_pos)=c("Monocytes","DCs Mig.","Macrophages","Neutrophils")
colnames(effect_non_spec)=c("Monocytes","DCs Mig.","NK","Lymphocytes","DC2s Res.","pDCs","DC1s Res.")


png("Project_Weizmann/NK_depletion_validation/Pagoda_2_results/Barplot_DE_genes_non_spec.png",
    width = 11,height = 12,units = "cm",res = 400,type = "cairo")
par(las=1,mar=c(7,4,2,2))
u=barplot(effect_non_spec,ylim=c(0,max((effect_non_spec))*1.3),col=c("black","white"),
          ylab="Number of DE genes",main="Non specific cells",beside=T,las=2,xlim=c(0,20))
dev.off()

png("Project_Weizmann/NK_depletion_validation/Pagoda_2_results/Barplot_DE_genes_ag_pos.png",
    width = 11,height = 12,units = "cm",res = 400,type = "cairo")
par(las=1,mar=c(7,4,2,2))
u=barplot(effect_ag_pos,ylim=c(0,max((effect_ag_pos))*1.3),col=c("black","white"),
          ylab="Number of DE genes",main="Antigen positive cells",beside=T,las=2,xlim=c(0,20))
dev.off()

#2)Volcanoplot 

mono_non_spec_logFC=aggregate(t(TPM_data_non_spec[,r_bis$clusters$PCA$infomap==1]),
                              by=list(factor(condition_non_spec_count[r_bis$clusters$PCA$infomap==1])),FUN=mean)
rownames(mono_non_spec_logFC)=mono_non_spec_logFC[,1]
mono_non_spec_logFC=t(mono_non_spec_logFC[,-1])
colnames(mono_non_spec_logFC)=c("Blocking","Control")
non_spec_loess=loess(Blocking~Control,data.frame(mono_non_spec_logFC),
                     degree = 2,span = 0.3,evaluation = 1000)

mono_non_spec_genes=rownames(DE_mono$table)[DE_mono$table$Lib_size>3]

par(las=1)

plot(non_spec_loess$residuals[mono_non_spec_genes],
     DE_mono$table[mono_non_spec_genes,"Condition"],
     xlim=c(-5,5),pch=16,
     col=string.to.colors(mono_non_spec_genes%in%c(rownames(DE_mono$table)[DE_mono$table$Condition>1.3]),
                          colors=c(alpha("red",0.7),alpha("black",0.3))),
     xlab="aNK1.1/Ifng vs Isotype (Log2 Fold Change)",ylab="-Log10(P-value)")
grid(nx = 4,ny=4,lty=2)



###
macro_infla_logFC=aggregate(t(TPM_data_ag_pos[,r$clusters$PCA$infomap==1]),
                            by=list(factor(condition_ag_pos_count[r$clusters$PCA$infomap==1])),FUN=mean)
rownames(macro_infla_logFC)=macro_infla_logFC[,1]
macro_infla_logFC=t(macro_infla_logFC[,-1])
colnames(macro_infla_logFC)=c("Blocking","Control")
macro_loess=loess(Blocking~Control,data.frame(macro_infla_logFC),
                  degree = 2,span = 0.3,evaluation = 1000)

macro_spec_genes=rownames(DE_macro_infla$table)[DE_macro_infla$table$Lib_size>3]

plot(macro_loess$residuals[macro_spec_genes],
     DE_macro_infla$table[macro_spec_genes,"Condition"],
     xlim=c(-6,6),pch=16,
     col=string.to.colors(macro_spec_genes%in%c(rownames(DE_macro_infla$table)[DE_macro_infla$table$Condition>1.3]),
                          col=(c(alpha("red",0.7),alpha("black",0.3)))),
     xlab="aNK1.1/Ifng vs Isotype (Log2 Fold Change)",ylab="-Log10(P-value)")
grid(nx = 4,ny=4,lty=2)


#VI)Reproductbility of the study

#A)Construction of the datasets

Mono_non_spec_reprod=aggregate(t(TPM_data_non_spec[,r_bis$clusters$PCA$infomap==1]),
                               by=list(replicates_non_spec_count[r_bis$clusters$PCA$infomap==1]),FUN=mean)
rownames(Mono_non_spec_reprod)=Mono_non_spec_reprod[,1]
Mono_non_spec_reprod=t(Mono_non_spec_reprod[,-1])
plot(Mono_non_spec_reprod[,1:2])

Mono_non_ag_pos_reprod=aggregate(t(TPM_data_ag_pos[,r$clusters$PCA$infomap==1]),
                                 by=list(replicates_ag_pos_count[r$clusters$PCA$infomap==1]),FUN=mean)
rownames(Mono_non_ag_pos_reprod)=Mono_non_ag_pos_reprod[,1]
Mono_non_ag_pos_reprod=t(Mono_non_ag_pos_reprod[,-1])
plot(Mono_non_ag_pos_reprod[,1:2])


#B)Plot of the results

par(las=1)
plot(Mono_non_spec_reprod[,1:2],pch=16,cex=1,
     xlab="Ms day2 + aNK1.1/Ifng rep.1",ylab="Ms day2 + aNK1.1/Ifng rep.2",
     xlim=c(0,16),ylim=c(0,16),col=alpha("black",0.2),main="Non specific Monocytes")
R=cor(Mono_non_spec_reprod[,1:2])[2,1]
legend("topleft",legend = paste("R = ",round(R,digits = 2),sep = ""),bty="n")
abline(lm(Mono_non_spec_reprod[,2]~Mono_non_spec_reprod[,1]),lwd=2,lty=2)

plot(Mono_non_spec_reprod[,3:4],pch=16,cex=1,
     xlab="Ms day2 + Isotype rep.1",ylab="Ms day2 + Isotype rep.2",
     xlim=c(0,16),ylim=c(0,16),col=alpha("black",0.2),main="Non specific Monocytes")
R=cor(Mono_non_spec_reprod[,3:4])[2,1]
legend("topleft",legend = paste("R = ",round(R,digits = 2),sep = ""),bty="n")
abline(lm(Mono_non_spec_reprod[,4]~Mono_non_spec_reprod[,3]),lwd=2,lty=2)

plot(Mono_non_ag_pos_reprod[,1:2],pch=16,cex=1,
     xlab="Ms day2 + aNK1.1/Ifng rep.1",ylab="Ms day2 + aNK1.1/Ifng rep.2",
     xlim=c(0,16),ylim=c(0,16),col=alpha("black",0.2),main="Antigen positive Monocytes")
R=cor(Mono_non_ag_pos_reprod[,1:2])[2,1]
legend("topleft",legend = paste("R = ",round(R,digits = 2),sep = ""),bty="n")
abline(lm(Mono_non_ag_pos_reprod[,2]~Mono_non_ag_pos_reprod[,1]),lwd=2,lty=2)

plot(Mono_non_ag_pos_reprod[,3:4],pch=16,cex=1,
     xlab="Ms day2 + Isotype rep.1",ylab="Ms day2 + Isotype rep.2",
     xlim=c(0,16),ylim=c(0,16),col=alpha("black",0.2),main="Antigen positive Monocytes")
R=cor(Mono_non_ag_pos_reprod[,3:4])[2,1]
legend("topleft",legend = paste("R = ",round(R,digits = 2),sep = ""),bty="n")
abline(lm(Mono_non_ag_pos_reprod[,4]~Mono_non_ag_pos_reprod[,3]),lwd=2,lty=2)





