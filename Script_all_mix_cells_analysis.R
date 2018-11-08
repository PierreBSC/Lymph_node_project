library(pagoda2)

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

#A)Non specific cells

data_final_1=read.table("Project_Weizmann/Final_sequencing/data_final.txt",row.names = 1,sep="\t")

annotation_1=read.table("Project_Weizmann/Final_sequencing/annotation.txt",row.names = 1,sep="\t")
annotation_1=as.character(annotation_1$x)
names(annotation_1)=colnames(data_final_1)

condition_1=read.table("Project_Weizmann/Final_sequencing/condition.txt",row.names = 1,sep="\t")
condition_1=as.character(condition_1$x)
names(condition_1)=colnames(data_final_1)

batch_1=read.table("Project_Weizmann/Final_sequencing/batch.txt",row.names = 1,sep="\t")
batch_1=as.character(batch_1$x)
names(batch_1)=colnames((data_final_1))


#B)Antigen specific cells

data_final_2=read.table("Project_Weizmann/Final_Antigen_pos/data_final.txt",row.names = 1,sep="\t")

condition_2=read.table("Project_Weizmann/Final_Antigen_pos/condition.txt",row.names = 1,sep="\t")
condition_2=as.character(condition_2$x)
names(condition_2)=colnames(data_final_2)

batch_2=read.table("Project_Weizmann/Final_Antigen_pos/batch.txt",row.names = 1,sep="\t")
batch_2=as.character(batch_2$x)
names(batch_2)=colnames(data_final_2)

#C)We merge the data

data_final=cbind(data_final_1,data_final_2)
condition=c(condition_1,condition_2)
batch=c(batch_1,batch_2)

antigen_status=c(rep("Non specific",length(batch_1)),rep("Antigen positive",length(batch_2)))
names(antigen_status)=colnames(data_final)

rm(data_final_1)
rm(data_final_2)


#####II)Filtering of the data

#A)Cell filtering
lib_size=colSums(data_final)
names(lib_size)=colnames(data_final)
par(las=1,family="serif")
hist(log10(1+lib_size),xlab="Library size of the cells (log10)",main="Librarys size distribution",n=100)
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

TPM_data=t(log2(1+t(as.matrix(data_count))/lib_size[colnames(data_count)]*10^6))
TPM_data=as.matrix(TPM_data)


#III)Analysis using Pagoda2 of the whole dataset

r <- Pagoda2$new(data_count,log.scale=FALSE)
r$adjustVariance(plot=T,gam.k=10)
r$calculatePcaReduction(nPcs=100,n.odgenes=3e3)
r$makeKnnGraph(k=40,type='PCA',center=T,distance='cosine')

r$getKnnClusters(method=multilevel.community,type='PCA')
r$getKnnClusters(method=infomap.community,type='PCA',name="infomap",)
r$getKnnClusters(method=walktrap.community,type='PCA',name="walktrap")


#### Two different clusters correspond to MigDCs clusters 8 and 22 and two to Monocytes (2 and 6) : we will focus on these cells for the later part of the analysis
##Of course the resultst names of the clustering will change according to the seed : change accordingly ...


DC_cells=names(r$clusters$PCA$community)[(r$clusters$PCA$community==8 |r$clusters$PCA$community==22)]
condition_DC=condition_count[DC_cells]
batch_DC=batch_count[DC_cells]
antigen_status_DC=antigen_status[DC_cells]

Mono_cells=names(r$clusters$PCA$community)[(r$clusters$PCA$community==2 |r$clusters$PCA$community==6)]
condition_mono=condition_count[Mono_cells]
batch_mono=batch_count[Mono_cells]
antigen_status_mono=antigen_status[Mono_cells]


#IV)DC analysis

#A)Analysis per se

r_DC <- Pagoda2$new(data_count[,DC_cells],log.scale=FALSE)
r_DC$adjustVariance(plot=T,gam.k=10)
r_DC$calculatePcaReduction(nPcs=100,n.odgenes=3e3)
r_DC$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')


r_DC$getKnnClusters(method=multilevel.community,type='PCA')
r_DC$getKnnClusters(method=infomap.community,type='PCA',name="infomap",)
r_DC$getKnnClusters(method=walktrap.community,type='PCA',name="walktrap")

r_DC$getEmbedding(type='PCA',embeddingType='largeVis',perplexity=20)
r_DC$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=20)


r_DC$plotEmbedding(type='PCA',embeddingType = "largeVis",groups = factor(antigen_status_DC),show.legend=F,mark.clusters=T,min.group.size=10,clusterType='community',
                   shuffle.colors=F,mark.cluster.cex=1,alpha=0.2,main='Phenograph clustering (T-SNE)')

r_DC$plotEmbedding(type='PCA',embeddingType = "tSNE",groups = factor(condition_DC),show.legend=F,mark.clusters=T,min.group.size=10,clusterType='community',
                   shuffle.colors=F,mark.cluster.cex=1,alpha=0.2,main='Phenograph clustering (T-SNE)')

#B)Over-dispersion analysis 

Gene_cor=cor(t(as.matrix(data_count[r_DC$getOdGenes(),names(which(selected_DC))])),method = "pearson")
Gene_clustering=hclust(dist(Gene_cor),method = "ward")
Gene_clustering=cutree(Gene_clustering,k = 20)
gene_env=c()
for (k in 1:length(unique(Gene_clustering))) {
  gene_env[[k]]=names(which(Gene_clustering==k))
}
names(gene_env)=1:20
gene_env_2 <- list2env(gene_env) # convert to an environment

r_DC$testPathwayOverdispersion(setenv = gene_env_2,type = "counts",verbose = T,
                               plot = F,max.pathway.size = 60,min.pathway.size = 1,recalculate.pca=T)

pathway_info=r_DC$misc$pathwayODInfo
pathway_info=pathway_info[order(pathway_info$cz,decreasing=T),]
barplot(pathway_info$cz)

dispersed_genes=c()

for (i in (substr(rownames(pathway_info),7,10))[1:20]) {
  u=r_DC$misc$pwpca[[i]]$xp$rotation
  u=(u[order(u[,1],decreasing = F),])
  dispersed_genes[[i]]=u
  print(i)
}
dispersed_genes=dispersed_genes[c(1,3,4,6:10)]

#C)Visualisation of the analysis 

PCA_score_DC=data.frame(MHC.II=as.numeric(r_DC$misc$pwpca$`15`$xp$scores),
                        Th2=as.numeric(r_DC$misc$pwpca$`3`$xp$scores),
                        Chimiokines=as.numeric(r_DC$misc$pwpca$`6`$xp$scores),
                        Th1=as.numeric(r_DC$misc$pwpca$`5`$xp$scores),
                        Cytoskeleton=as.numeric(r_DC$misc$pwpca$`10`$xp$scores),
                        Costim_1=as.numeric(r_DC$misc$pwpca$`17`$xp$scores),
                        Costim_2=as.numeric(r_DC$misc$pwpca$`7`$xp$scores),
                        Migration=as.numeric(r_DC$misc$pwpca$`14`$xp$scores),
                        row.names =colnames(r_DC$misc$pwpca$`14`$xp$scores) )
PCA_score_DC=PCA_score_DC[selected_DC,]

names_pathways_DC=c("MHC-II genes","DC activation","Chimiokines","Th1 cytokines",
                    "Cytoskeleton","Costimulation pathway","Costimulation pathway 2","Migration")

color_ordered=col=alpha(c("grey","salmon","salmon","darkred","darkred",
                          "cornflowerblue","cornflowerblue","darkblue","darkblue",
                          "palegreen","palegreen","darkgreen","darkgreen"),alpha = 0.7)
names_ordered=c("PBS","Ms day1 Ag+","Ms day1 Ag-","Ms day2 Ag+","Ms day2 Ag-",
                "Nb day1 Ag+","Nb day1 Ag-","Nb day2 Ag+","Nb day2 Ag-",
                "Ca day1 Ag+","Ca day1 Ag-","Ca day2 Ag+","Ca day2 Ag-")


mixed_condition_DC_bis=paste(condition_DC_bis,antigen_status_DC_bis,sep = "")
l=unique(mixed_condition_DC_bis)[order(unique(mixed_condition_DC_bis))]
mixed_condition_DC_bis=factor(mixed_condition_DC_bis,levels = l[c(5:13,1:4)])

pdf("Project_Weizmann/Final_all_mix/Pagoda_2_analysis/DC_score.pdf",width = 10,height = 6)

for (k in 1:length(dispersed_genes)) {
  split.screen(rbind(c(0,0.3,0,1),c(0.3,1,0,1)))
  screen(1)
  par(las=1,mar=c(6,5,2,2))
  barplot(dispersed_genes[[k]],horiz = T,
          col="black",xlab="Gene contribution to the PCA",
          xlim=c(0,max(dispersed_genes[[k]])*1.2))
  screen(2)
  par(las=2,mar=c(8,5,2,2))
  boxplot(PCA_score_DC[,k]~mixed_condition_DC_bis,outline=F,
          col=color_ordered,names=names_ordered,ylab="PCA Score",
          main=names_pathways_DC[k],cex.lab=1.5)
  close.screen(all.screens = TRUE)
}
dev.off()


#V)Mono analysis

#A)Analysis per se 

r_mono <- Pagoda2$new(data_count[,Mono_cells],log.scale=FALSE)
r_mono$adjustVariance(plot=T,gam.k=10)
r_mono$calculatePcaReduction(nPcs=100,n.odgenes=3e3)
r_mono$makeKnnGraph(k=30,type='PCA',center=T,distance='cosine')

r_mono$getKnnClusters(method=multilevel.community,type='PCA')
r_mono$getKnnClusters(method=infomap.community,type='PCA',name="infomap")
r_mono$getKnnClusters(method=walktrap.community,type='PCA',name="walktrap")

#B)Over-dispersion analysis 

Gene_cor=cor(t(as.matrix(data_count[r_mono$getOdGenes(),names(which(selected_mono))])),method = "pearson")
Gene_clustering=hclust(dist(Gene_cor),method = "ward")
Gene_clustering=cutree(Gene_clustering,k = 30)
gene_env=c()
for (k in 1:length(unique(Gene_clustering))) {
  gene_env[[k]]=names(which(Gene_clustering==k))
}
names(gene_env)=1:30
gene_env_2 <- list2env(gene_env) # convert to an environment

r_mono$testPathwayOverdispersion(setenv = gene_env_2,type = "counts",verbose = T,
                                 plot = F,max.pathway.size = 60,min.pathway.size = 1,recalculate.pca=T)

pathway_info=r_mono$misc$pathwayODInfo
pathway_info=pathway_info[order(pathway_info$cz,decreasing=T),]
barplot(pathway_info$cz)

dispersed_genes=c()

for (i in (substr(rownames(pathway_info),7,10))[1:15]) {
  u=r_mono$misc$pwpca[[i]]$xp$rotation
  u=u[unique(rownames(u)),]
  u=(u[order(u,decreasing = F)])
  dispersed_genes[[i]]=u
  print(i)
}
dispersed_genes=dispersed_genes[c(1,2,4,6,8,9,11)]

#C)Visualisation of the analysis 

PCA_score_mono=data.frame(Ifnb=as.numeric(r_mono$misc$pwpca$`27`$xp$scores),
                          Ifng=as.numeric(r_mono$misc$pwpca$`5`$xp$scores),
                          Cathepsin=as.numeric(r_mono$misc$pwpca$`8`$xp$scores),
                          MHC.II=as.numeric(r_mono$misc$pwpca$`20`$xp$scores),
                          M2=as.numeric(r_mono$misc$pwpca$`13`$xp$scores),
                          Th1=as.numeric(r_mono$misc$pwpca$`17`$xp$scores),
                          C1=as.numeric(r_mono$misc$pwpca$`15`$xp$scores),
                          row.names=colnames(r_mono$misc$pwpca$`5`$xp$scores))
PCA_score_mono=PCA_score_mono[selected_mono,]

mixed_condition_mono_bis=paste(condition_mono_bis,antigen_status_mono_bis,sep = "")
l=unique(mixed_condition_mono_bis)[order(unique(mixed_condition_mono_bis))]
mixed_condition_mono_bis=factor(mixed_condition_mono_bis,levels = l[c(5:13,1:4)])

names_pathways_mono=c(expression(paste("IFN",beta," pathway",sep = "")),
                      expression(paste("LPS/IFN",gamma," pathway",sep = "")),
                      "Cathepsin pathway","MHC-II genes",
                      "M2 polarisation","Th1 cytokines","C1q genes")

color_ordered=col=alpha(c("grey","salmon","salmon","darkred","darkred",
                          "cornflowerblue","cornflowerblue","darkblue","darkblue",
                          "palegreen","palegreen","darkgreen","darkgreen"),alpha = 0.7)
names_ordered=c("PBS","Ms day1 Ag+","Ms day1 Ag-","Ms day2 Ag+","Ms day2 Ag-",
                "Nb day1 Ag+","Nb day1 Ag-","Nb day2 Ag+","Nb day2 Ag-",
                "Ca day1 Ag+","Ca day1 Ag-","Ca day2 Ag+","Ca day2 Ag-")

pdf("Project_Weizmann/Final_all_mix/Pagoda_2_analysis/Mono_score.pdf",width = 10,height = 6)

for (k in 1:length(dispersed_genes)) {
  split.screen(rbind(c(0,0.3,0,1),c(0.3,1,0,1)))
  screen(1)
  par(las=1,mar=c(6,5,2,2))
  barplot(dispersed_genes[[k]],horiz = T,
          col="black",xlab="Gene contribution to the PCA",
          xlim=c(0,max(dispersed_genes[[k]])*1.2))
  screen(2)
  par(las=2,mar=c(7,5,2,2))
  boxplot(PCA_score_mono[,k]~mixed_condition_mono_bis,outline=F,
          col=color_ordered,names=names_ordered,ylab="PCA Score",
          main=names_pathways_mono[k],cex.lab=1.5)
  close.screen(all.screens = TRUE)
}
dev.off()



##VI)Violin plots for monocytes and DCs across conditions

condition_DC_bis=factor(condition_DC_bis,levels = c("Control","Ms day1","Ms day2",
                                                    "Nb day1","Nb day2",
                                                    "Ca day1",'Ca day2'))
antigen_status_DC_bis=factor(antigen_status_DC_bis,levels = c("Antigen positive","Non specific"))

#A)For DC

violin_plot_DC=function(gene) {
  data_gene=data.frame(Expression=TPM_data_DC[gene,],
                       Condition=condition_DC_bis,
                       Ag_status=antigen_status_DC_bis)
  ggplot(data_gene, aes(x=Condition, y=Expression,fill=Condition,linetype=Ag_status))  +  geom_violin(trim=T,scale = "width",bw=1.5,draw_quantiles = 0.5,na.rm = T,show.legend = F)+
    scale_y_continuous(name = "Expression log2(TPM)",limits = c(0,max(TPM_data_DC[gene,]))) + scale_fill_manual(values = color_ordered_bis)+
    theme_classic() + ggtitle(gene) + theme(plot.title = element_text(size=22),axis.text=element_text(size = 15),axis.title = element_text(size = 15)) + scale_x_discrete(labels=NULL,name=" ")  
}
par(las=1,mfrow=c(4,2))

pdf("Project_Weizmann/Revision_plot/Violin_plot_all_mix_DC.pdf",width = 5,height = 3.5)
violin_plot_DC("Ccl22")
violin_plot_DC("Ccl17")
violin_plot_DC("Cd40")
violin_plot_DC("Cd86")
dev.off()

pdf("Project_Weizmann/Revision_plot/Violin_plot_all_mix_DC_bis.pdf",width = 5,height = 3.5)
violin_plot_DC("Ccl3")
violin_plot_DC("Cxcl2")
violin_plot_DC("Cxcl3")
dev.off()

#B)For Monocytes

condition_mono_bis=factor(condition_mono_bis,levels = c("Control","Ms day1","Ms day2",
                                                        "Nb day1","Nb day2",
                                                        "Ca day1",'Ca day2'))
antigen_status_mono_bis=factor(antigen_status_mono_bis,levels = c("Antigen positive","Non specific"))

color_ordered_bis=unique(color_ordered)

violin_plot_mono=function(gene) {
  data_gene=data.frame(Expression=TPM_data_mono[gene,],
                       Condition=condition_mono_bis,
                       Ag_status=antigen_status_mono_bis)
  ggplot(data_gene, aes(x=Condition, y=Expression,fill=Condition,linetype=Ag_status))  +  geom_violin(trim=T,scale = "width",bw=1.5,draw_quantiles = 0.5,na.rm = T,show.legend = F)+
    scale_y_continuous(name = "Expression log2(TPM)",limits = c(0,max(TPM_data_mono[gene,]))) + scale_fill_manual(values = color_ordered_bis)+
    theme_classic() + ggtitle(gene) + theme(plot.title = element_text(size=22),axis.text=element_text(size = 15),axis.title = element_text(size = 15)) + scale_x_discrete(labels=NULL,name=" ")  
}

pdf("Project_Weizmann/Revision_plot/Violin_plot_all_mix_mono.pdf",width = 5,height = 3.5)
violin_plot_mono("Il12b")
violin_plot_mono("Tnf")
violin_plot_mono("Ctsd")
violin_plot_mono("Mrc1")
dev.off()

pdf("Project_Weizmann/Revision_plot/Violin_plot_all_mix_mono_sup.pdf",width = 5,height = 3.5)
violin_plot_mono("Cd40")
violin_plot_mono("C3ar1")
violin_plot_mono("Ctsb")
dev.off()
