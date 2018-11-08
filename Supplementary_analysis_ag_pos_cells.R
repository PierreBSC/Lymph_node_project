
### Warning : this script can only be used if the script for the analysis ofAntigen positive cells has been run before



##I)Detailled QC of the data

#A)Library size


#General distribution
pdf("Project_Weizmann/Revision_plot/QC_Ag_pos.pdf",width = 5.5,height = 5.5)
lib_size=colSums(data_final)
names(lib_size)=colnames(data_final)
par(las=1)
hist(log10(lib_size),xlab="Library size of the cells (Log10)",main="Library size distribution",
     n=100,xlim=c(1,5),xaxs='i',yaxs='i')
abline(v=log10(350),lwd=2,lty=2)

#Across sequencing batches
batch=batch[colnames(data_final)]
boxplot(log10(lib_size)~batch,outline=F,xlab="Sequencing batches",
        ylab="Library size of the cells (Log10)",cex.axis=0.7,cex.lab=1.3)
abline(h=log10(350),lwd=2,lty=2)

#Across clusters
lib_size_clusters=lib_size[colnames(heatmap_data)]
boxplot(log10(lib_size_clusters)~factor(r$clusters$PCA$community[colnames(heatmap_data)],
                                        levels = order_cluster_bis),outline=F,names=1:6,
        ylab="Library size of the cells (Log10)",xlab="Cluster",
        main="Library size across clusters",cex.lab=1.3)



#Across samples

cell_proportion=table(replicates_count,r$clusters$PCA$community)
cell_proportion=cell_proportion[c(6:13,1,2,4,5),order_cluster_bis]
sample_names=rownames(cell_proportion)
lib_size_replicates=lib_size[replicates%in%sample_names]
selected_samples=factor(replicates[replicates%in%sample_names],sample_names)


color_ordered=col=alpha(c("salmon","salmon","darkred","darkred",
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
boxplot(log10(Number_gene_detected)~batch,outline=F,xlab="Sequencing batches",
        ylab="Number of detected genes (Log10)")

plot(log10(lib_size),log10(Number_gene_detected),
     xlab="Library size (Log10)",ylab="Number of genes (Log10)",
     pch=16,col=alpha("black",0.2),main="Relation between library size and number of genes detected")
abline(v=log10(350),lwd=2,lty=2)


##Across cluster

Number_gene_detected_cluster=Number_gene_detected[colnames(heatmap_data)]
boxplot(log10(Number_gene_detected_cluster)~factor(r$clusters$PCA$community[colnames(heatmap_data)],
                                                   levels = order_cluster_bis),outline=F,names=1:6,
        ylab="Number of genes (Log10)",xlab="Cluster",
        main="Number of genes detected  across clusters")

#C)Proportion of ERCC genes

ERCC_genes=rownames(data_final)
ERCC_genes=ERCC_genes[grepl(pattern = "ERCC",ERCC_genes)]
ERCC_UMI=colSums(data_final[ERCC_genes,])
ERCC_UMI=ERCC_UMI/lib_size

hist(ERCC_UMI,xlab="Percent of ERCC among total UMIs",n=100,main="")
boxplot(ERCC_UMI~batch,outline=F,xlab="Sequencing batches",
        ylab="Percent of ERCC among total UMIs")

ERCC_cluster=ERCC_UMI[colnames(heatmap_data)]
boxplot((ERCC_cluster)~factor(r$clusters$PCA$community[colnames(heatmap_data)],
                              levels = order_cluster_bis),outline=F,names=1:6,
        ylab="Percent of ERCC among total UMIs",xlab="Cluster",
        main="Quality of the cells across clusters",ylim=c(0,1))

#D)Contribution of each sample for the analysis

v=table(replicates_count)[c(6:13,1,2,4,5)]
barplot(v,xlab="Replicates",ylab="Number of filtered cells",col=color_ordered)

dev.off()

##III)FACS data : can we predict DC transcriptional status solely based on Index Sorting data

FACS_data=read.table("Project_Weizmann/Final_Antigen_pos/FACS_data_final.txt")
FACS_data=FACS_data[colnames(data_count),]
rownames(FACS_data)=colnames(data_count)

plot(FACS_data$FSC.A,FACS_data$FceRIa,pch=21,bg="red",cex=0.5)

DC_cells=r$clusters$PCA$community%in%c(8,9)
FACS_data_DC=FACS_data[DC_cells,]
FACS_data_DC=FACS_data_DC[,c(-1,-16,-17)]
DC_cluster=r$clusters$PCA$community[DC_cells]
FACS_data_DC$Cluster=ifelse(DC_cluster==unique(DC_cluster)[1],1,0)


FACS_data_DC=FACS_data_DC[,c(1,4,7,9:11,13:15)]

DC_model=glm(Cluster~.,FACS_data_DC,family = "binomial")
summary(DC_model)
get_ROC_curve=function(logistic_fit,true_values,main="ROC curve SSA") {
  prob=predict(logistic_fit,type=c("response"))    
  library("ROCR")    
  pred <- prediction(prob, true_values)    
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")     
  plot(perf, col=rainbow(7), main=main, xlab="Specificity", 
       ylab="Sensitivity",lwd=2.5,cex.lab=1.5,cex.main=1.5,cex.axis=1.5)    
  abline(0, 1,lwd=2.5) #add a 45 degree line
  perf <- performance(pred, measure = "auc")
  legend("bottomright",bty="n",cex=1.5,
         legend = paste("AUC =",round(unlist(perf@y.values),digits = 2)))
}


coef_logis_classifier=summary(DC_model)
coef_logis_classifier=coef_logis_classifier$coefficients[-1,4]
coef_logis_classifier=-log10(coef_logis_classifier)
coef_logis_classifier=coef_logis_classifier[order(coef_logis_classifier,decreasing = F)]

get_ROC_curve(logistic_fit = DC_model,true_values = FACS_data_DC$Cluster[complete.cases(FACS_data_DC)],
              main = "Prediction of DC subtype")
barplot(coef_logis_classifier,horiz = T,xlim=c(0,2.5),xlab="-Log10(P-value)",
        main="Variable contribution to logistic classifier ")



