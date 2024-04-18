#tutorial from https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html
#install.packages("vegan")
library("tidyverse")
library("vegan")
library(dplyr)
library(stringr)
library(data.table)        
#install.packages("FactoMineR")
library("FactoMineR")
#install.packages("factoextra")
library(factoextra)
#BiocManager::install("edgeR")
library(edgeR)
library(limma)
#BiocManager::install("snm")
library(snm)
#BiocManager::install("Maaslin2")
library(Maaslin2)
#BiocManager::install("M3C")
library(M3C)
library("pvca")
library(lme4)
library("labdsv")
library("coin")
library("vegan")
library("yaml")
library("ggpubr")
library("cowplot")
library("tidyverse")
library("dplyr")
library("UpSetR")


#Upload metadata
metadata=read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples_QC/GTEx_Argis_QCed_all_metadata.tab")

#Upload agamenon results data
#Read agamemnon results files of all the GTEx 
GTEx_all_run = read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples/ml_models_after_bowtie/species/ML_model/agamemnon_after_bowtie_allSamples_species.tab")
colnames(GTEx_all_run)=gsub("\\.", "-",colnames(GTEx_all_run))
metadata$samples.submitter_id
which(colnames(GTEx_all_run)%in%metadata$samples.submitter_id)
all=GTEx_all_run[,which(colnames(GTEx_all_run)%in%c("External-ID",metadata$samples.submitter_id))]    
rownames(GTEx_all_run)

#Remove the first column of all and make the rownames of all the different species
rownames=all$"External-ID"
all=all[,2:length(colnames(all))] #Columns the samples and rows the microorganisms
rownames(all)=rownames
dim(all)
colnames(all)

#From all remove the species with zero value in all the samples 
all_sub=all[rowSums(all==0, na.rm=TRUE)<ncol(all), ]

#Change the order of the records of metadata according to all_sub columns
metadata2=metadata[which(metadata$samples.submitter_id%in%colnames(all_sub)), ]
metadata2=metadata2[match(colnames(all_sub), metadata2$samples.submitter_id), ]
dim(metadata2)

#alpha diversity 
Richiness=colSums(all_sub != 0)
invsimpson_diversity= diversity(all_sub,
                                MARGIN = 2,
                                index= "invsimpson")
simpson_diversity= diversity(all_sub,
                             MARGIN = 2,
                             index= "simpson")
shannon_diversity= diversity(all_sub,
                             MARGIN = 2,
                             index= "shannon")

data_diversity=data.frame(Richiness,invsimpson_diversity, simpson_diversity,shannon_diversity,tissue_type=metadata2$tissue_type, sampleID=metadata2$samples.submitter_id)

# anova for shannon diveristy
anova=aov(shannon_diversity~tissue_type, data=data_diversity)
summary(anova)

#Post hoc analysis
tukey <- TukeyHSD(anova)
print(tukey)

# anova for richineess
anova=aov(Richiness~tissue_type, data=data_diversity)
summary(anova)

#Post hoc analysis
tukey <- TukeyHSD(anova)
print(tukey)


#Creating plots
data_diversity2=data_diversity %>% mutate(tissue_type = tissue_type) %>%
  group_by(tissue_type) %>% 
  mutate(mean_shannon=mean(shannon_diversity), mean_richiness=mean(Richiness),
         mean_simpson=mean(simpson_diversity), mean_invSimpson=mean(invsimpson_diversity))

jpeg(file="alpha_diverisity_per_tissue_all.jpeg",width=20,height=10,units="in",res=600)
ggplot(data_diversity2,aes(x=tissue_type, y=shannon_diversity, fill=mean_shannon)) +
  geom_boxplot(fatten = NULL) +
  theme_classic() + 
  stat_summary(fun.y=mean, geom="point", color="black", aes(group=tissue_type)) +
  ylab('Alpha diversity (Shannon)') +
  xlab('') +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 34),
        axis.text = element_text(size = 28),
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position= 'none')
dev.off()

jpeg(file="richiness_per_tissue_all_github.jpeg",width=25,height=10,units="in",res=600)
ggplot(data_diversity2,aes(x=tissue_type, y=Richiness, fill=mean_richiness)) +
  geom_boxplot(fatten = NULL) +
  theme_bw() + 
  stat_summary(fun.y=mean, geom="point", color="black", aes(group=tissue_type)) +
  ylab('Richness') +
  xlab('') + 
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        text = element_text(size = 34),
        axis.text = element_text(size = 28),
        axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position= 'none')  
dev.off()


#################################################################################################################################################
################################################ Find core microbiome with Filter10Percent #######################################################################################
#Microbiomes with at least 1 read at at least 10 percent of samples of at least one tissue

#Core microbiome function
find_core_microbiome_Filter10= function(data_frame, metadata, tissue_type) {
  tissue_samples=metadata[which(metadata$tissue_type==tissue_type),"specimen_id"]
  data_test=as.data.frame(data_frame)
  data_tissue=data_test[,tissue_samples]
  core_data_tissue=data_tissue[rowSums(data_tissue>1, na.rm=TRUE)>(0.1*ncol(data_tissue)), ]
  core_data_tissue_microbiomes=rownames(core_data_tissue)
  return(core_data_tissue_microbiomes)
}

#Create a variable named core_"tissue" per tissue, in which there is saved the core microbiome of this tissue
core=c()
core_dataframe=data.frame(matrix(ncol=1))
for (tissue in unique(metadata2$tissue_type)) {
  print(tissue)
  core_microbiome_name=paste("core_",tolower(tissue),sep="")
  core_microbiome_name=gsub(" ", "_", core_microbiome_name)
  print(core_microbiome_name)
  assign(core_microbiome_name, find_core_microbiome_Filter10(all_sub,metadata2,tissue))
  core=c(core,find_core_microbiome_Filter10(all_sub,metadata2,tissue))
  all_sub_tissue=all_sub[,which(metadata2$tissue_type==tissue)]
  all_sub_tissue$mean_frequency=rowMeans(all_sub_tissue)
  core_tissue=find_core_microbiome_Filter10(all_sub,metadata2,tissue)
  core_dataframe=qpcR:::cbind.na(core_dataframe, core_tissue,
                                 all_sub_tissue$mean_frequency[which(rownames(all_sub_tissue)%in%core_tissue)])
}

core=unique(core)
core_dataframe=core_dataframe[,2:length(colnames(core_dataframe))]
colnames_used=c()
for (name in unique(metadata2$tissue_type)) {
  colnames_used=c(colnames_used,name,paste(name,"_meanFrequency",sep=""))
} 
colnames(core_dataframe)=colnames_used

# "Methanosarcina sp. WH1"%in%core
# core_dataframe_muscle=core_dataframe[,c("Muscle", "Muscle_meanFrequency")]
# core_dataframe_muscle[which(core_dataframe_muscle$Muscle=="Methanosarcina sp. WH1"),]
# core_dataframe_bladder=core_dataframe[,c("Bladder", "Bladder_meanFrequency")]
# core_dataframe_bladder[which(core_dataframe_bladder$Bladder=="Caldimonas thermodepolymerans"),]
# core_dataframe_bl=core_dataframe[,c("Bladder", "Bladder_meanFrequency")]
# core_dataframe_bladder[which(core_dataframe_bladder$Bladder=="Caldimonas thermodepolymerans"),]
# core_dataframe_bl=core_dataframe[,c("Blood", "Blood_meanFrequency")]
# core_dataframe_bl[which(core_dataframe_bl$Blood=="Caldimonas thermodepolymerans"),]

write.csv(core_dataframe, file="core_microbiome_per_tissue.csv", row.names = F)


#Find the components for each kigndom in the core microbiome
kingdom_species=read.delim("kingdom_species.txt",sep=",",header=FALSE)
bacteria=kingdom_species[which(kingdom_species$V1=="Bacteria"),]
archaea=kingdom_species[which(kingdom_species$V1=="Archaea"),]
fungi=kingdom_species[which(kingdom_species$V1=="Eukaryota"),]
virus=kingdom_species[-which(kingdom_species$V2%in%fungi$V2 | kingdom_species$V2%in%bacteria$V2), ]

#Find the core only for bacteria 
core_bacteria=c()
core_dataframe_bacteria=data.frame(matrix(ncol=1))
for (tissue in unique(metadata2$tissue_type)) {
  print(tissue)
  core_microbiome_name=paste("core_",tolower(tissue),sep="")
  core_microbiome_name=gsub(" ", "_", core_microbiome_name)
  print(core_microbiome_name)
  assign(core_microbiome_name, find_core_microbiome_Filter10(all_sub,metadata2,tissue))
  core_temp = find_core_microbiome_Filter10(all_sub,metadata2,tissue)
  core_temp=core_temp[core_temp%in%bacteria$V2]
  core_dataframe_bacteria=qpcR:::cbind.na(core_dataframe_bacteria, core_temp)
  core_bacteria=c(core_bacteria,core_temp)
}
core_bacteria=unique(core_bacteria)
core_dataframe_bacteria=core_dataframe_bacteria[,2:length(colnames(core_dataframe_bacteria))]
colnames(core_dataframe_bacteria)=unique(metadata2$tissue_type)
write.csv(core_dataframe_bacteria, file="core_microbiome_per_tissue_bacteria.csv", row.names = F)



#Find the core only for fungi 
core_fungi=c()
core_dataframe_fungi=data.frame(matrix(ncol=1))
for (tissue in unique(metadata2$tissue_type)) {
  print(tissue)
  core_microbiome_name=paste("core_",tolower(tissue),sep="")
  core_microbiome_name=gsub(" ", "_", core_microbiome_name)
  print(core_microbiome_name)
  assign(core_microbiome_name, find_core_microbiome_Filter10(all_sub,metadata2,tissue))
  core_temp = find_core_microbiome_Filter10(all_sub,metadata2,tissue)
  core_temp=core_temp[core_temp%in%fungi$V2]
  core_dataframe_fungi=qpcR:::cbind.na(core_dataframe_fungi, core_temp)
  core_fungi=c(core_fungi,core_temp)
}
core_fungi=unique(core_fungi)
core_dataframe_fungi=core_dataframe_fungi[,2:length(colnames(core_dataframe_fungi))]
colnames(core_dataframe_fungi)=unique(metadata2$tissue_type)
write.csv(core_dataframe_fungi, file="core_microbiome_per_tissue_fungi.csv", row.names = F)



#Find the core only for Virus 
core_virus=c()
core_dataframe_virus=data.frame(matrix(ncol=1))
for (tissue in unique(metadata2$tissue_type)) {
  print(tissue)
  core_microbiome_name=paste("core_",tolower(tissue),sep="")
  core_microbiome_name=gsub(" ", "_", core_microbiome_name)
  print(core_microbiome_name)
  assign(core_microbiome_name, find_core_microbiome_Filter10(all_sub,metadata2,tissue))
  core_temp = find_core_microbiome_Filter10(all_sub,metadata2,tissue)
  core_temp=core_temp[core_temp%in%virus$V2]
  core_dataframe_virus=qpcR:::cbind.na(core_dataframe_virus, core_temp)
  core_virus=c(core_virus,core_temp)
}
core_virus=unique(core_virus)
core_dataframe_virus=core_dataframe_virus[,2:length(colnames(core_dataframe_virus))]
colnames(core_dataframe_virus)=unique(metadata2$tissue_type)
write.csv(core_dataframe_virus, file="core_microbiome_per_tissue_virus.csv", row.names = F)


#Find the core only for archae 
core_archae=c()
core_dataframe_archae=data.frame(matrix(ncol=1))
for (tissue in unique(metadata2$tissue_type)) {
  print(tissue)
  core_microbiome_name=paste("core_",tolower(tissue),sep="")
  core_microbiome_name=gsub(" ", "_", core_microbiome_name)
  print(core_microbiome_name)
  assign(core_microbiome_name, find_core_microbiome_Filter10(all_sub,metadata2,tissue))
  core_temp = find_core_microbiome_Filter10(all_sub,metadata2,tissue)
  core_temp=core_temp[core_temp%in%archaea$V2]
  core_dataframe_archae=qpcR:::cbind.na(core_dataframe_archae, core_temp)
  core_archae=c(core_archae,core_temp)
}
core_archae=unique(core_archae)
core_dataframe_archae=core_dataframe_archae[,2:length(colnames(core_dataframe_archae))]
colnames(core_dataframe_archae)=unique(metadata2$tissue_type)
write.csv(core_dataframe_archae, file="core_microbiome_per_tissue_archaea.csv", row.names = F)

#Create a dataframe in order to create the barplots of kingdom distrbution of each tissue's core microbiome
dataframe_kingdoms_core=data.frame(matrix(nrow=4, ncol=28))
colnames(dataframe_kingdoms_core)= unique(metadata2$tissue_type)
rownames(dataframe_kingdoms_core)=c("Bacteria","Archaea","Fungi","Viruses")
for (tissue in unique(metadata2$tissue_type)){
  name=paste("core_",tolower(tissue),sep="")
  name=gsub(" ", "_", name)
  core_temp=get(paste0(name))
  dataframe_kingdoms_core[,tissue]=c(length(which(core_temp%in%bacteria$V2)),length(which(core_temp%in%archaea$V2)),
                                     length(which(core_temp%in%fungi$V2)),length(which(core_temp%in%virus$V2)))
}

#Make the plots 
# Define the number of times to repeat each element
repeat_factor <- 4
# Repeat each element in the first list and combine them into a new list
ColumnName_list <- rep(unique(metadata2$tissue_type), each = repeat_factor)

# Convert original dataframe to a new dataframe with two columns
new_dataframe_kingdoms_core <- data.frame(
  Value = unlist(dataframe_kingdoms_core),
  ColumnName = ColumnName_list,
  kingdom = unlist(rep(c("Bacteria","Archaea","Fungi","Viruses"),28))
)

#Make Histograms
count=seq(1,28*5,step=1)
count=count[unlist(count) %% 5 != 0]
length(count)
dim(new_dataframe_kingdoms_core)
p1=ggplot(transform(new_dataframe_kingdoms_core, x=count), aes(x = ColumnName, y = Value, fill = kingdom)) +
  geom_bar(stat = "identity", position="dodge") +
  labs(title = "Histogram of species' kingdoms per tissue", x = "", y = "Number of species") +
  geom_text(aes(label=Value), position = position_dodge(width = 1), vjust = -0.2, size=4) +
  theme_bw() +
  scale_fill_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size=20), text = element_text(size = 24)) 
p2=ggplot(new_dataframe_kingdoms_core, aes(x = ColumnName, y = Value, fill = kingdom)) +
  geom_bar(stat = "identity", position="dodge") +
  labs(x = "", y = "Number of species") +
  geom_text(aes(label=Value), position = position_dodge(width = 1), vjust = -0.2, size=6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size=25), text = element_text(size = 25)) 


jpeg(file="histogram_perKingdom_perTissue2_core.jpg", width=20, height=10, units="in", res=450)
p2
dev.off()



#Upset diagrams for the shared total core microbiome 
my_list_core=as.list(core_dataframe)
my_list_core = rapply(my_list_core, na.omit, how = "replace")

jpeg(file="Upset_core_microbiome.jpeg",width=12,height=5,units="in",res=600)
upset(fromList(my_list_core), nsets=30, order.by = "degree",
      mainbar.y.label = "Tissue Intersections", sets.x.label = "Core Microbiome",
      mb.ratio = c(0.45, 0.55))
dev.off()


#Upset diagrams for the shared bacterial core microbiome 
my_list_bacteria=as.list(core_dataframe_bacteria)
my_list_bacteria = rapply(my_list_bacteria, na.omit, how = "replace")


jpeg(file="Upset_core_microbiome_bacteria.jpeg",width=12,height=5,units="in",res=600)
upset(fromList(my_list_bacteria), nsets=30, order.by = "degree",
      mainbar.y.label = "Tissue Intersections", sets.x.label = "Core Bacterial Microbiome",
      mb.ratio = c(0.45, 0.55))
dev.off()

#Upset diagrams for the shared fungal core microbiome 
my_list_fungi=as.list(core_dataframe_fungi)
my_list_fungi = rapply(my_list_fungi, na.omit, how = "replace")


jpeg(file="Upset_core_microbiome_fungi.jpeg",width=12,height=5,units="in",res=600)
upset(fromList(my_list_fungi), nsets=30, order.by = "degree",
      mainbar.y.label = "Tissue Intersections", sets.x.label = "Core Fungal Microbiome",
      mb.ratio = c(0.45, 0.55))
dev.off()


#Upset diagrams for the shared viral core microbiome 
my_list_virus=as.list(core_dataframe_virus)
my_list_virus = rapply(my_list_virus, na.omit, how = "replace")

jpeg(file="Upset_core_microbiome_virus.jpeg",width=12,height=5,units="in",res=600)
upset(fromList(my_list_virus), nsets=30, order.by = "degree",
      mainbar.y.label = "Tissue Intersections", sets.x.label = "Core Viral Microbiome",
      mb.ratio = c(0.45, 0.55))
dev.off()

#Save in csv files shared core bactreria, fungi and virus
shared_bacteria <- Reduce(intersect, my_list_bacteria)
write.table(shared_bacteria, file="shared_core_microbiome_bacteria_allTissues.csv", col.names=FALSE, sep=",")

shared_fungi <- Reduce(intersect, my_list_fungi)
write.table(shared_fungi, file="shared_core_microbiome_fungi_allTissues.csv", col.names=FALSE, sep=",")

shared_viral <- Reduce(intersect, my_list_virus)
write.table(shared_viral, file="shared_core_microbiome_virus_allTissues.csv", col.names=FALSE, sep=",")

