`%nin%`=Negate(`%in%`)
library(rsample)      # data splitting
# model visualization
library(dplyr)
library(ggplot2)      # model visualization
# 
#Run the commented code one time in order to create the data.RData and then run the script from the command line 
# 
# #Upload metadata
# metadata=read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples_QC/GTEx_Argis_QCed_all_metadata.tab")
# dim(metadata)
# #Upload agamenon results data
# 
# #Read agamemnon results file for the demo run
# GTEx_all_run = read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples/ml_models_after_bowtie/species/ML_model/agamemnon_after_bowtie_allSamples_species.tab")
# colnames(GTEx_all_run)=gsub("\\.", "-",colnames(GTEx_all_run))
# metadata$samples.submitter_id
# which(colnames(GTEx_all_run)%in%metadata$samples.submitter_id)
# all=GTEx_all_run[,which(colnames(GTEx_all_run)%in%c("External-ID",metadata$samples.submitter_id))]    
# 
# #Remove the first column of all and make the rownames of all the different species
# rownames=all$"External-ID"
# all=all[,2:length(colnames(all))] #Columns the samples and rows the microorganisms
# rownames(all)=rownames
# dim(all)
# colnames(all)
# 
# #From all remove the species with zero value in all the samples 
# all_sub=all[rowSums(all==0, na.rm=TRUE)<ncol(all), ]
# 
# 
# #Change the order of the records of metadata according to all_sub columns
# metadata2=metadata[which(metadata$samples.submitter_id%in%colnames(all_sub)), ]
# metadata2=metadata2[match(colnames(all_sub), metadata2$samples.submitter_id), ]
# rownames(metadata2)=metadata2$samples.submitter_id
# dim(metadata2)
# 
# #################################################################################################################################################
# ################################################ Normalise the datasets with CSS  #######################################################################################
# 
# library(metagenomeSeq)
# metaSeqObject = newMRexperiment(all_sub) 
# metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
# dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
# residuals=dge_normalised
# residuals[is.na(residuals)]=0
# colnames(residuals)=colnames(all_sub)
# 
# save(residuals,all_sub,metadata2, file = "data.RData")

load("data.RData")
#Input the microbiome wanted
args = commandArgs(trailingOnly = TRUE)
microbiome=args[1]
normalised=args[2]
#microbiome="Allomeiothermus silvanus"
#normalised="CSS"
print(microbiome)
print(normalised)

#data=t(all_sub)
#data=t(residuals)
if (normalised == "CSS"){
  data=t(residuals)
} else if (normalised == "Raw"){
  data=t(all_sub)
} else {
  print("Error in normised argument")
}
# rownames(data)
# rownames(metadata2)

#Create a dataframe for the boxplot
if (length(which(rownames(data)==rownames(metadata2)))==length(rownames(data))){
  dataframe_microbiome=data.frame("microbiome"=data[,microbiome], "subjects"=metadata2$samples.submitter_id, "tissue_type"=metadata2$tissue_type)
}
#Plot

dataframe_microbiome2=dataframe_microbiome %>%  group_by(tissue_type) %>% mutate(mean_microbiome = mean(microbiome))

jpeg(file=paste(microbiome,"_tissueExpression_",normalised,".jpeg",sep=""),width=15,height=8,units="in",res=150)

ggplot(dataframe_microbiome2,aes(x=tissue_type, y=microbiome, fill=mean_microbiome)) +
  geom_boxplot(fatten = NULL) +
  theme_classic() + 
  stat_summary(fun.y=mean, geom="point", color="black", aes(group=tissue_type)) +
  ylab(paste('Normalised Counts of\n',microbiome,sep="")) +
  xlab('') +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position= 'none')  
dev.off()
