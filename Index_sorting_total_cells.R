#I) Selection of the cells
source("~/AFA_pipeline_script.R") # This house script contains all the scripts needed for INX/Flow cytometry extraction

#A)Loading and selection 

annotation=read.delim("Project_Weizmann/UMI.tab/Correspondance.txt",header=TRUE,sep="\t")

#First step "Selection of the interesting samples : 
#migth be either B/T or non B/T but no stromal cells and no Antigen cells 

annotation_bis=annotation[grep(pattern = "non B/T",x = annotation$gating..cells.sorted.),]

# We select only the experiments with sequenced data
annotation_bis=annotation_bis[annotation_bis$Batch..!="",]

# We remove the VSV and CMV experiments plus the "none" experiment 
annotation_bis=annotation_bis[-c(grep(pattern = c("CMV"),annotation_bis$treatment.of.mice),
                                 grep(pattern = c("VSV"),annotation_bis$treatment.of.mice),
                                 grep(pattern = c("Isotype"),annotation_bis$treatment.of.mice,ignore.case = T),
                                 grep(pattern = c("anti"),annotation_bis$treatment.of.mice)),]


# We select only the experiments with sequenced data
annotation_bis=annotation_bis[annotation_bis$Batch..!="",]

list_plaques=as.character(annotation_bis$Plate..)
list_plaques=paste(list_plaques,".fcs",sep = "")
kept_plaques=list_plaques%in%(list.files("Project_Weizmann/FCS_files/")) ## Which plaques do we really have....
list_plaques=list_plaques[kept_plaques]

specific_files=which((list.files("Project_Weizmann/FCS_files/"))%in%list_plaques)

###II)Loading, scaling, normalisation of the data using the home-made pipeline

non_spec_data=AFA_pipe_load(loading_dir = "Project_Weizmann/FCS_files/",specific_file = specific_files,change_name = F)
non_spec_data=AFA_pipe_transform(non_spec_data,chanel_to_normalise = 8:16)


for (k in 1:length(non_spec_data)) {
  comp <- keyword(non_spec_data[[k]])$`SPILL`
  non_spec_data[[k]] <- compensate(non_spec_data[[k]],comp)
  
}
correspondance_channel=pData(parameters(non_spec_data[[1]]))
new_chanel_names=ifelse(is.na(correspondance_channel$desc),yes = correspondance_channel$name,no = correspondance_channel$desc)

colnames(non_spec_data)=as.character(new_chanel_names)
colnames(non_spec_data)=make.names(colnames(non_spec_data))


###III)Extraction of the Index sorting data 

#A)Creation of a clean correspondance file

correspondance=read.table("Project_Weizmann/wells_cells.txt",header=TRUE,sep="\t",row.names = 1)
correspondance=correspondance[correspondance$Amp_batch_ID%in%annotation_bis$Batch..,]

correspondance=data.frame(Cell_ID=rownames(correspondance),
                          Position=correspondance$well_coordinates,
                          Batch=correspondance$Amp_batch_ID,
                          Plate=rep(NA,length(correspondance$Amp_batch_ID)))


for (k in unique(as.character(correspondance$Batch))) {
  u=annotation_bis[annotation_bis$Batch..==k,1]
  correspondance[correspondance$Batch==k,"Plate"]=u
}
correspondance$Plate=paste(correspondance$Plate,".fcs",sep = "")

#B)Extraction of the index sorting data

INX_data=AFA_index_sorting(flow_data = non_spec_data,correspondance)
rownames(INX_data)=INX_data$Cell_ID
INX_data=INX_data[,-ncol(INX_data)]

