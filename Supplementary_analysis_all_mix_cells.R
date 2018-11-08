### Warning : this script can only be used if the script for the analysis of all mixed  cells has been run before

##I)Study of DCs in details 

##DCs :cluster 8 and 22

r$getDifferentialGenes(type = "PCA",clusterType = "infomap",z.threshold = 3,upregulated.only = T,verbose = T)

count_DC_subtype=table(as.character(r$clusters$PCA$infomap[r$clusters$PCA$community%in%c(8,22)]))
count_DC_subtype=count_DC_subtype[order(count_DC_subtype,decreasing = T)]
selected_DC_pop=names(count_DC_subtype[1:10])
count_DC_subtype=table(as.character(r$clusters$PCA$infomap[r$clusters$PCA$infomap%in%selected_DC_pop]))
count_DC_subtype=count_DC_subtype[selected_DC_pop]

gene_list=c()

for (k in selected_DC_pop) {
  u=r$diffgenes$PCA$infomap[[k]] 
  u=u[u$highest,]
  u=u[order(u$M,decreasing = T),]
  gene_list=c(gene_list,rownames(u[1:5,]))
}
gene_list=gene_list[!grepl("NA",gene_list)]

heatmap_DC=c()

for (k in as.numeric(selected_DC_pop)) {
  u=TPM_data[gene_list,r$clusters$PCA$infomap==k]
  heatmap_DC=cbind(heatmap_DC,u)
}

heatmap_DC=as.matrix(heatmap_DC)
lim_values=quantile(as.numeric((heatmap_DC)),c(0.05,0.99))

heatmap_DC[heatmap_DC<lim_values[1]]=lim_values[1]
heatmap_DC[heatmap_DC>lim_values[2]]=lim_values[2]

png("Project_Weizmann/Revision_plot/DC_heatmap.png", width=6000, height=6000,type = "cairo")
par(las=1,family="Arial",mar=c(130,40,4,4))
image(rotate(heatmap_DC[,]),xaxt="n",yaxt='n',col=colorRampPalette(c("white","white","white","gold","orange","maroon4"))(100))
abline(v=cumsum(count_DC_subtype[1:10])/sum(count_DC_subtype[1:10]),lwd=15,lty=2)
box(which = "plot",lty = 1,lwd=20)
axis(2, at=seq(0,1,length.out=nrow((heatmap_DC[,]))), 
     labels= ( rownames(heatmap_DC)[nrow(heatmap_DC):1]),
     las= 2, cex.axis=7, tick = F )
dev.off()
##
DC_list=which(r$clusters$PCA$infomap%in%selected_DC_pop)
DC_enrichment=table(condition_count[DC_list],as.character(r$clusters$PCA$infomap[DC_list]))

DC_enrichment=round(DC_enrichment/rowSums(DC_enrichment)*100,digits = 0)
DC_enrichment=DC_enrichment[,selected_DC_pop]
DC_enrichment=DC_enrichment["Control",]

pdf("Project_Weizmann/Revision_plot/barplot_DC_control",width = 8,height = 4)
barplot(matrix(DC_enrichment),horiz = T,
        xlab="Proportion of cluster among control samples",col = rainbow(10))
dev.off()

get_DC_expression=function(gene) {
  u=r$counts[r$clusters$PCA$infomap%in%selected_DC_pop,gene]
  v=factor(r$clusters$PCA$infomap[r$clusters$PCA$infomap%in%selected_DC_pop],selected_DC_pop)
  boxplot(u~v)
}


