`%nin%`=Negate(`%in%`)
library(rsample)      # data splitting 
library(gbm)          # basic implementation
library(xgboost)      # a faster implementation of gbm
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # a java-based platform
library(pdp)          # model visualization
library(ggplot2)      # model visualization
library(lime)         # model visualization
# Loading package
library(dplyr)
# Loading package
library(caTools)
library("pROC")
library(ROCR)
#install.packages("PRROC")
library("PRROC")
#Libraries Needed
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
#PVCA install
#https://rdrr.io/github/dleelab/pvca/#vignettes
#remotes::install_github("dleelab/pvca", force=TRUE)
library("pvca")
library(lme4)
#install_github("vqv/ggbiplot")
#install_github("road2stat/ggsci")
library(devtools)
library(ggsci)
#if (!requireNamespace("remotes", quietly=TRUE))
#  install.packages("remotes")
#remotes::install_github("yiluheihei/microbiomeMarker")
library(microbiomeMarker)
library(metagenomeSeq)
library(wesanderson)


#############################################################################################################################
######################################### HEALTHY DATASET ###################################################################

#Upload healhy dataset
data_healthy_first=read.delim("/mnt/raid1/argis/GTEx/after_access/normal_samples_OtherDatesets/tissues/PRJEB4337/Agamemnon_after_bowtie/final/all/Agamemnon_after_bowtie_species_PRJEB44337.tab")
metadata_healthy=readRDS("/mnt/raid1/argis/GTEx/after_access/normal_samples_OtherDatesets/metadata_healthy.Rdata")
rownames(data_healthy_first)=data_healthy_first$External.ID
data_healthy_first=data_healthy_first[,2:length(colnames(data_healthy_first))]
#Normalise with CSS
metaSeqObject = newMRexperiment(data_healthy_first) 
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
data_normalised_healthy_first=dge_normalised
data_normalised_healthy_first[is.na(data_normalised_healthy_first)]=0

#Check halthy metadata
metadata_healthy=metadata_healthy[which(metadata_healthy$sample%in%colnames(data_normalised_healthy_first)),]
rownames(metadata_healthy)=seq(1, length(rownames(metadata_healthy)), by=1)

#metadata_healthy$project=c(rep("Random",219),rep("Project",171))
table(metadata_healthy$tissue)
table(metadata_healthy[,2:3])

#Keep only the samples for the 8 tissues 
metadata_healthy=metadata_healthy[which(tolower(metadata_healthy$tissue)%in%c("heart","liver", "bladder", "muscle",
                                                                              "stomach","colon","testis", "blood")),]
data_normalised_healthy_first=data_normalised_healthy_first[,which(colnames(data_normalised_healthy_first)%in%metadata_healthy$sample)]

dim(data_normalised_healthy_first)
dim(metadata_healthy)
table(metadata_healthy$tissue)
data_normalised_healthy_first[1:5,1:5]

# Convert values to relative frequency
data_relative_healthy = as.data.frame(lapply(data_normalised_healthy_first, function(x) x / sum(x)))
data_relative_tissues=data.frame(matrix(nrow=length(data_normalised_healthy_first[,1]), ncol =1))
for  (tissue in unique(metadata_healthy$tissue)) {
  data_relative_health_tissue=data_relative_healthy[,which(metadata_healthy$tissue%in%tissue)]
  data_relative_tissues=cbind(data_relative_tissues,as.data.frame( rowMeans(data_relative_health_tissue)))
}
data_relative_tissues=data_relative_tissues[,2:length(colnames(data_relative_tissues))]
colnames(data_relative_tissues)=unique(metadata_healthy$tissue)
rownames(data_relative_tissues)=rownames(data_normalised_healthy_first)
dim(data_relative_tissues)

write.csv(data_relative_tissues,file="Relative_Frequency_per_tissue_living_dataset.csv")

# For liver
data_relative_tissues$microbiomes= rownames(data_relative_tissues)
liver_relative=as.data.frame(data_relative_tissues[,c("Liver","microbiomes")],)
liver_relative_sorted = as.data.frame(liver_relative[order(liver_relative[,1],decreasing = TRUE), ])
features = rownames(liver_relative_sorted)[1:5]
values = liver_relative_sorted[1:5,1]  # Example values for the features

# Create a dataframe with features and their corresponding values
data_liver = data.frame(features, values)
data_liver$features = str_wrap(data_liver$features, width = 2)  # Adjust width as needed

pie_chart = ggplot(data_liver, aes(x = "", y = values, fill = features)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(x = 1.3, label = paste(features, "\n", round(values*100),"%"), angle=0, size =14), 
            position = position_stack(vjust = 0.5)) +
  coord_polar("y") +
  labs(title = "", fill = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "none") + 
  scale_fill_manual(values=c(wes_palette("GrandBudapest2")[2], wes_palette("Darjeeling2")[4], 
                             wes_palette("Darjeeling2")[2],wes_palette("FantasticFox1")[2],
                             wes_palette("Darjeeling1")[4]))
pie_chart
# Show the plot

# Create pie chart with ggplot2
png("liver_healthy_relative_frequnecy.png", width = 7.5, height = 7.5, units = "in", res=600)
pie_chart
dev.off()    



# For Colon
data_relative_tissues$microbiomes= rownames(data_relative_tissues)
colon_relative=as.data.frame(data_relative_tissues[,c("Colon","microbiomes")],)
colon_relative_sorted = as.data.frame(colon_relative[order(colon_relative[,1],decreasing = TRUE), ])
features = rownames(colon_relative_sorted)[1:5]
values = colon_relative_sorted[1:5,1]  # Example values for the features

# Create a dataframe with features and their corresponding values
data_colon = data.frame(features, values)
data_colon$features = str_wrap(data_colon$features, width = 2)  # Adjust width as needed

pie_chart = ggplot(data_colon, aes(x = "", y = values, fill = features)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(x = 1.3, label = paste(features, "\n", round(values*100),"%"), angle=0, size =14), 
            position = position_stack(vjust = 0.5)) +
  coord_polar("y") +
  labs(title = "", fill = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "none")+ 
  scale_fill_manual(values=c(wes_palette("GrandBudapest2")[2],wes_palette("Darjeeling2")[4],  wes_palette("Darjeeling2")[2], 
                             wes_palette("Moonrise2")[3], wes_palette("Darjeeling1")[4]))
# Show the plot

# Create pie chart with ggplot2
png("colon_healthy_relative_frequnecy.png", width = 7.5, height = 7.5, units = "in", res=600)
pie_chart
dev.off()    



# For Heart
data_relative_tissues$microbiomes= rownames(data_relative_tissues)
heart_relative=as.data.frame(data_relative_tissues[,c("Heart","microbiomes")],)
heart_relative_sorted = as.data.frame(heart_relative[order(heart_relative[,1],decreasing = TRUE), ])
features = rownames(heart_relative_sorted)[1:5]
values = heart_relative_sorted[1:5,1]  # Example values for the features

# Create a dataframe with features and their corresponding values
data_heart = data.frame(features, values)
data_heart$features = c(str_wrap(data_heart$features[1], width = 2), data_heart$features[2:5])  # Adjust width as needed

pie_chart = ggplot(data_heart, aes(x = "", y = values, fill = features)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(x = 1.3, label = paste(features,"(", round(values*100),"%)"), angle=50, size =14), 
            position = position_stack(vjust = 0.5)) +
  coord_polar("y") +
  labs(title = "", fill = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "none")+ 
  scale_fill_manual(values=c(wes_palette("GrandBudapest2")[2],wes_palette("Darjeeling2")[4],  wes_palette("Darjeeling1")[4],
                             wes_palette("Darjeeling2")[2], wes_palette("Moonrise2")[3]))

pie_chart
# Show the plot

# Create pie chart with ggplot2
png("heart_healthy_relative_frequnecy.png", width = 7.5, height = 7.5, units = "in", res=600)
pie_chart
dev.off()    
