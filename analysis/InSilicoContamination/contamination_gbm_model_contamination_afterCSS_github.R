library(rsample)      # data splitting 
library(gbm)          # basic implementation
library(xgboost)      # a faster implementation of gbm
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # a java-based platform
library(pdp)          # model visualization
#library(ggplot2)      # model visualization
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
#library("pvca")
library(lme4)
#install_github("vqv/ggbiplot")
#install_github("road2stat/ggsci")
library(devtools)
library(ggsci)
if (!requireNamespace("remotes", quietly=TRUE))
  install.packages("remotes")
#remotes::install_github("yiluheihei/microbiomeMarker")
library(microbiomeMarker)
library(metagenomeSeq)


#Upload metadata
metadata=read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples_QC/GTEx_Argis_QCed_all_metadata.tab")
dim(metadata)

#Upload agamenon results data
#Read agamemnon results file for the demo run
GTEx_all_run = read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples/ml_models_after_bowtie/species/ML_model/agamemnon_after_bowtie_allSamples_species.tab")
colnames(GTEx_all_run)=gsub("\\.", "-",colnames(GTEx_all_run))
all=GTEx_all_run[,which(colnames(GTEx_all_run)%in%c("External-ID",metadata$samples.submitter_id))]    

#Remove the first column of all and make the rownames of all the different species
rownames=all$"External-ID"
all=all[,2:length(colnames(all))] #Columns the samples and rows the microorganisms
rownames(all)=rownames

#From all remove the species with zero value in all the samples 
all_sub=all[rowSums(all==0, na.rm=TRUE)<ncol(all), ]


#Change the order of the records of metadata according to all_sub columns
metadata2=metadata[which(metadata$samples.submitter_id%in%colnames(all_sub)), ]
metadata2=metadata2[match(colnames(all_sub), metadata2$samples.submitter_id), ]
dim(metadata2)


#Keep the 11 tissues with tissue-specific microbiome
tissues_keep=c("heart","liver","small intestine", "bladder", "brain", "muscle",
               "stomach","colon","testis", "blood","salivary gland")
all_sub=all_sub[,which(tolower(metadata2$tissue_type)%in%c("heart","liver","small intestine", "bladder", "brain", "muscle",
                                                           "stomach","colon","testis", "blood","salivary gland")),]

dim(all_sub)
dim(metadata2)
metadata2=metadata2[which(tolower(metadata2$tissue_type)%in%c("heart","liver","small intestine", "bladder", "brain", "muscle",
                                                              "stomach","colon","testis", "blood","salivary gland")),]
table(metadata2$tissue_type)


#################################################################################################################################################
################################################ Find core microbiome with Filter10Percent #######################################################################################

#Microbiomes with filter that for each microbiome there is at least 1 read at at least 10 percent of samples of at least one tissue

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
for (tissue in unique(metadata2$tissue_type)) {
  print(tissue)
  core_microbiome_name=paste("core_",tolower(tissue),sep="")
  core_microbiome_name=gsub(" ", "_", core_microbiome_name)
  print(core_microbiome_name)
  assign(core_microbiome_name, find_core_microbiome_Filter10(all_sub,metadata2,tissue))
  core=c(core,find_core_microbiome_Filter10(all_sub,metadata2,tissue))
}

core=unique(core)


####################################################################################
################### For each time of the 100 itterations ###########################
set.seed(42)
for (i in seq(1,100,step=1)) {
  #i=5
  
  ####################################################################################################################################
  ################################## Split the datasets to train and testing dataset 
  index <- createDataPartition(metadata2$tissue_type, p = 0.7, list = FALSE)
  trainX <- all_sub[,index]
  metadata_train <- metadata2[index,]
  testX <- all_sub[,-index]
  metadata_test <- metadata2[-index,]
  dim(trainX)
  dim(testX)
  table(metadata_train$tissue_type)
  table(metadata_test$tissue_type)
  dim(all_sub)
  
  #################################################################################################################################################
  ################################################ Separately normalise the datasets with CSS Maaslin2 #######################################################################################
  
  #Train dataset
  metaSeqObject = newMRexperiment(trainX) 
  metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
  dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
  residuals_train=dge_normalised
  residuals_train[is.na(residuals_train)]=0
  colnames(residuals_train)=colnames(trainX)
  residuals_train2=residuals_train[core,]
  
  
  #Train dataset
  metaSeqObject = newMRexperiment(testX) 
  metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
  dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
  residuals_test=dge_normalised
  residuals_test[is.na(residuals_test)]=0
  rownames(residuals_test)
  colnames(residuals_test)=colnames(testX)
  residuals_test2=residuals_test[core,]
  
  contamination_file_txt=paste("./contamination/test_",i,"_contamination.csv",sep="")
  output_file=paste("./results/test_",i,"_permutations_results.log",sep="")
  contamination_file=read.csv(contamination_file_txt)
  
  ###################################################################################################################################
  ################################# Create Contaminataion values #####################################################################
  #Create contamination for residuals_test2
  all_contamin=residuals_test2[1:2, ]
  #Insert the contaminations 
  contamination_1=rep(0,length(colnames(residuals_test2)))
  #contamination_1_train=rep(0,length(colnames(residuals_train2)))
  all_1 <- rbind(all_contamin, contamination_1)
  rownames(all_1)=c(rownames(all_contamin), "contamination_1")
  # Remove the surrounding square brackets
  contamination_1_list_string = gsub("\\[|\\]", "", contamination_file[1,2])
  # Remove the surrounding single quotes
  string_contamination_1 <- gsub("'", "", contamination_1_list_string)
  # Split the string into a list of strings
  string_list_contamination_1 <- strsplit(string_contamination_1, ", ")[[1]]
  all_1[which(rownames(all_1)=="contamination_1"), which(colnames(all_1)%in%string_list_contamination_1)]=150
  
  
  contamination_2=rep(0,length(colnames(all_1)))
  all_2 <- rbind(all_1, contamination_2)
  rownames(all_2)=c(rownames(all_1), "contamination_2")
  contamination_2_list_string = gsub("\\[|\\]", "", contamination_file[2,2])
  string_contamination_2 <- gsub("'", "", contamination_2_list_string)
  string_list_contamination_2 <- strsplit(string_contamination_2, ", ")[[1]]
  all_2[which(rownames(all_2)=="contamination_2"), which(colnames(all_2)%in%string_list_contamination_2)]=150
  
  contamination_3=rep(0,length(colnames(all_2)))
  all_3 <- rbind(all_2, contamination_3)
  rownames(all_3)=c(rownames(all_2), "contamination_3")
  contamination_3_list_string = gsub("\\[|\\]", "", contamination_file[3,2])
  string_contamination_3 <- gsub("'", "", contamination_3_list_string)
  string_list_contamination_3 <- strsplit(string_contamination_3, ", ")[[1]]
  all_3[which(rownames(all_3)=="contamination_3"), which(colnames(all_3)%in%string_list_contamination_3)]=150
  
  contamination_4=rep(0,length(colnames(all_3)))
  all_4 <- rbind(all_3, contamination_4)
  rownames(all_4)=c(rownames(all_3), "contamination_4")
  contamination_4_list_string = gsub("\\[|\\]", "", contamination_file[4,2])
  string_contamination_4 <- gsub("'", "", contamination_4_list_string)
  string_list_contamination_4 <- strsplit(string_contamination_4, ", ")[[1]]
  all_4[which(rownames(all_4)=="contamination_4"), which(colnames(all_4)%in%string_list_contamination_4)]=150
  
  contamination_5=rep(0,length(colnames(all_4)))
  all_5 <- rbind(all_4, contamination_5)
  rownames(all_5)=c(rownames(all_4), "contamination_5")
  contamination_5_list_string = gsub("\\[|\\]", "", contamination_file[5,2])
  string_contamination_5 <- gsub("'", "", contamination_5_list_string)
  string_list_contamination_5 <- strsplit(string_contamination_5, ", ")[[1]]
  all_5[which(rownames(all_5)=="contamination_5"), which(colnames(all_5)%in%string_list_contamination_5)]=150
  
  
  contamination_6=rep(0,length(colnames(all_5)))
  all_6 <- rbind(all_5, contamination_6)
  rownames(all_6)=c(rownames(all_5), "contamination_6")
  contamination_6_list_string = gsub("\\[|\\]", "", contamination_file[6,2])
  string_contamination_6 <- gsub("'", "", contamination_6_list_string)
  string_list_contamination_6 <- strsplit(string_contamination_6, ", ")[[1]]
  all_6[which(rownames(all_6)=="contamination_6"), which(colnames(all_6)%in%string_list_contamination_6)]=150
  
  contamination_file[7,2]
  contamination_7=rep(0,length(colnames(all_6)))
  all_7 <- rbind(all_6, contamination_7)
  rownames(all_7)=c(rownames(all_6), "contamination_7")
  contamination_7_list_string = gsub("\\[|\\]", "", contamination_file[7,2])
  string_contamination_7 <- gsub("'", "", contamination_7_list_string)
  string_list_contamination_7 <- strsplit(string_contamination_7, ", ")[[1]]
  all_7[which(rownames(all_7)=="contamination_7"), which(colnames(all_7)%in%string_list_contamination_7)]=1
  
  
  contamination_8=rep(0,length(colnames(all_7)))
  all_8 <- rbind(all_7, contamination_8)
  rownames(all_8)=c(rownames(all_7), "contamination_8")
  contamination_8_list_string = gsub("\\[|\\]", "", contamination_file[8,2])
  string_contamination_8 <- gsub("'", "", contamination_8_list_string)
  string_list_contamination_8 <- strsplit(string_contamination_8, ", ")[[1]]
  #string_list_contamination_8=metadata[which(metadata$bss_collection_site=="C1"), ]$samples.submitter_id
  all_8[which(rownames(all_8)=="contamination_8"), which(colnames(all_8)%in%string_list_contamination_8)]=1
  
  contamination_9=rep(0,length(colnames(all_8)))
  all_9 <- rbind(all_8, contamination_9)
  rownames(all_9)=c(rownames(all_8), "contamination_9")
  contamination_9_list_string = gsub("\\[|\\]", "", contamination_file[9,2])
  string_contamination_9 <- gsub("'", "", contamination_9_list_string)
  string_list_contamination_9 <- strsplit(string_contamination_9, ", ")[[1]]
  all_9[which(rownames(all_9)=="contamination_9"), which(colnames(all_9)%in%string_list_contamination_9)]=1
  
  contamination_10=rep(0,length(colnames(all_9)))
  all_10 <- rbind(all_9, contamination_10)
  rownames(all_10)=c(rownames(all_9), "contamination_10")
  contamination_10_list_string = gsub("\\[|\\]", "", contamination_file[10,2])
  string_contamination_10 <- gsub("'", "", contamination_10_list_string)
  string_list_contamination_10 <- strsplit(string_contamination_10, ", ")[[1]]
  all_10[which(rownames(all_10)=="contamination_10"), which(colnames(all_10)%in%string_list_contamination_10)]=1
  
  contamination_11=rep(0,length(colnames(all_10)))
  all_11 <- rbind(all_10, contamination_11)
  rownames(all_11)=c(rownames(all_10), "contamination_11")
  contamination_11_list_string = gsub("\\[|\\]", "", contamination_file[11,2])
  string_contamination_11 <- gsub("'", "", contamination_11_list_string)
  string_list_contamination_11 <- strsplit(string_contamination_11, ", ")[[1]]
  all_11[which(rownames(all_11)=="contamination_11"), which(colnames(all_11)%in%string_list_contamination_11)]=1
  
  contamination_12=rep(0,length(colnames(all_11)))
  all_12 <- rbind(all_11, contamination_12)
  rownames(all_12)=c(rownames(all_11), "contamination_12")
  contamination_12_list_string = gsub("\\[|\\]", "", contamination_file[12,2])
  string_contamination_12 <- gsub("'", "", contamination_12_list_string)
  string_list_contamination_12 <- strsplit(string_contamination_12, ", ")[[1]]
  all_12[which(rownames(all_12)=="contamination_12"), which(colnames(all_12)%in%string_list_contamination_12)]=1
  contaminants_test2=tail(all_12,12)
  dim(contaminants_test2)
  residuals_test3=rbind(residuals_test2,contaminants_test2)
  
  
  #Create contamination for residuals_train2
  all_contamin=residuals_train2[1:2, ]
  #Insert the contaminations 
  contamination_1=rep(0,length(colnames(residuals_train2)))
  #contamination_1_train=rep(0,length(colnames(residuals_train2)))
  all_1 <- rbind(all_contamin, contamination_1)
  rownames(all_1)=c(rownames(all_contamin), "contamination_1")
  # Remove the surrounding square brackets
  contamination_1_list_string = gsub("\\[|\\]", "", contamination_file[1,2])
  # Remove the surrounding single quotes
  string_contamination_1 <- gsub("'", "", contamination_1_list_string)
  # Split the string into a list of strings
  string_list_contamination_1 <- strsplit(string_contamination_1, ", ")[[1]]
  all_1[which(rownames(all_1)=="contamination_1"), which(colnames(all_1)%in%string_list_contamination_1)]=150
  
  
  
  contamination_2=rep(0,length(colnames(all_1)))
  all_2 <- rbind(all_1, contamination_2)
  rownames(all_2)=c(rownames(all_1), "contamination_2")
  contamination_2_list_string = gsub("\\[|\\]", "", contamination_file[2,2])
  string_contamination_2 <- gsub("'", "", contamination_2_list_string)
  string_list_contamination_2 <- strsplit(string_contamination_2, ", ")[[1]]
  all_2[which(rownames(all_2)=="contamination_2"), which(colnames(all_2)%in%string_list_contamination_2)]=150
  #all_2[which(rownames(all_2)=="contamination_2"), 31:60]
  #"GTEX-11EQ9-1826-SM-5Q5AJ"%in%string_list_contamination_2
  
  contamination_3=rep(0,length(colnames(all_2)))
  all_3 <- rbind(all_2, contamination_3)
  rownames(all_3)=c(rownames(all_2), "contamination_3")
  contamination_3_list_string = gsub("\\[|\\]", "", contamination_file[3,2])
  string_contamination_3 <- gsub("'", "", contamination_3_list_string)
  string_list_contamination_3 <- strsplit(string_contamination_3, ", ")[[1]]
  all_3[which(rownames(all_3)=="contamination_3"), which(colnames(all_3)%in%string_list_contamination_3)]=150
  
  contamination_4=rep(0,length(colnames(all_3)))
  all_4 <- rbind(all_3, contamination_4)
  rownames(all_4)=c(rownames(all_3), "contamination_4")
  contamination_4_list_string = gsub("\\[|\\]", "", contamination_file[4,2])
  string_contamination_4 <- gsub("'", "", contamination_4_list_string)
  string_list_contamination_4 <- strsplit(string_contamination_4, ", ")[[1]]
  all_4[which(rownames(all_4)=="contamination_4"), which(colnames(all_4)%in%string_list_contamination_4)]=150
  
  contamination_5=rep(0,length(colnames(all_4)))
  all_5 <- rbind(all_4, contamination_5)
  rownames(all_5)=c(rownames(all_4), "contamination_5")
  contamination_5_list_string = gsub("\\[|\\]", "", contamination_file[5,2])
  string_contamination_5 <- gsub("'", "", contamination_5_list_string)
  string_list_contamination_5 <- strsplit(string_contamination_5, ", ")[[1]]
  all_5[which(rownames(all_5)=="contamination_5"), which(colnames(all_5)%in%string_list_contamination_5)]=150
  
  
  contamination_6=rep(0,length(colnames(all_5)))
  all_6 <- rbind(all_5, contamination_6)
  rownames(all_6)=c(rownames(all_5), "contamination_6")
  contamination_6_list_string = gsub("\\[|\\]", "", contamination_file[6,2])
  string_contamination_6 <- gsub("'", "", contamination_6_list_string)
  string_list_contamination_6 <- strsplit(string_contamination_6, ", ")[[1]]
  all_6[which(rownames(all_6)=="contamination_6"), which(colnames(all_6)%in%string_list_contamination_6)]=150
  
  contamination_file[7,2]
  contamination_7=rep(0,length(colnames(all_6)))
  all_7 <- rbind(all_6, contamination_7)
  rownames(all_7)=c(rownames(all_6), "contamination_7")
  contamination_7_list_string = gsub("\\[|\\]", "", contamination_file[7,2])
  string_contamination_7 <- gsub("'", "", contamination_7_list_string)
  string_list_contamination_7 <- strsplit(string_contamination_7, ", ")[[1]]
  all_7[which(rownames(all_7)=="contamination_7"), which(colnames(all_7)%in%string_list_contamination_7)]=1
  
  
  contamination_8=rep(0,length(colnames(all_7)))
  all_8 <- rbind(all_7, contamination_8)
  rownames(all_8)=c(rownames(all_7), "contamination_8")
  contamination_8_list_string = gsub("\\[|\\]", "", contamination_file[8,2])
  string_contamination_8 <- gsub("'", "", contamination_8_list_string)
  string_list_contamination_8 <- strsplit(string_contamination_8, ", ")[[1]]
  all_8[which(rownames(all_8)=="contamination_8"), which(colnames(all_8)%in%string_list_contamination_8)]=1
  
  contamination_9=rep(0,length(colnames(all_8)))
  all_9 <- rbind(all_8, contamination_9)
  rownames(all_9)=c(rownames(all_8), "contamination_9")
  contamination_9_list_string = gsub("\\[|\\]", "", contamination_file[9,2])
  string_contamination_9 <- gsub("'", "", contamination_9_list_string)
  string_list_contamination_9 <- strsplit(string_contamination_9, ", ")[[1]]
  all_9[which(rownames(all_9)=="contamination_9"), which(colnames(all_9)%in%string_list_contamination_9)]=1
  
  contamination_10=rep(0,length(colnames(all_9)))
  all_10 <- rbind(all_9, contamination_10)
  rownames(all_10)=c(rownames(all_9), "contamination_10")
  contamination_10_list_string = gsub("\\[|\\]", "", contamination_file[10,2])
  string_contamination_10 <- gsub("'", "", contamination_10_list_string)
  string_list_contamination_10 <- strsplit(string_contamination_10, ", ")[[1]]
  all_10[which(rownames(all_10)=="contamination_10"), which(colnames(all_10)%in%string_list_contamination_10)]=1
  
  contamination_11=rep(0,length(colnames(all_10)))
  all_11 <- rbind(all_10, contamination_11)
  rownames(all_11)=c(rownames(all_10), "contamination_11")
  contamination_11_list_string = gsub("\\[|\\]", "", contamination_file[11,2])
  string_contamination_11 <- gsub("'", "", contamination_11_list_string)
  string_list_contamination_11 <- strsplit(string_contamination_11, ", ")[[1]]
  all_11[which(rownames(all_11)=="contamination_11"), which(colnames(all_11)%in%string_list_contamination_11)]=1
  
  contamination_12=rep(0,length(colnames(all_11)))
  all_12 <- rbind(all_11, contamination_12)
  rownames(all_12)=c(rownames(all_11), "contamination_12")
  contamination_12_list_string = gsub("\\[|\\]", "", contamination_file[12,2])
  string_contamination_12 <- gsub("'", "", contamination_12_list_string)
  string_list_contamination_12 <- strsplit(string_contamination_12, ", ")[[1]]
  all_12[which(rownames(all_12)=="contamination_12"), which(colnames(all_12)%in%string_list_contamination_12)]=1
  contaminants_train2=tail(all_12,12)
  dim(contaminants_train2)
  residuals_train3=rbind(residuals_train2,contaminants_train2)
  dim(residuals_train3)
  
  #############################################################################################################################
  ################################ MACHINE LEARNING MODELS ####################################################################
  #############################################################################################################################
  #Parallel computing
  library(parallel)
  library(doMC) # for parallel computing
  numCores <- detectCores()
  #Set how many cores the script will use (10 cores)
  registerDoMC(cores=8)
  
  
  #The dataframe that will keep all the results of the tissues that have aML with good performance 
  results_dataframe=data.frame(matrix(nrow=1,ncol=length(tissues_keep)*14))
  column_names=c()
  for (tissue in unique(metadata2$tissue_type)) {
    column_names=c(column_names,paste(tissue,"_ROC",sep=""),paste(tissue,"_PR",sep=""))
    for (contaminant in rownames(contaminants_test2)) {
      column_names=c(column_names,paste(tissue,"_FIscore_",contaminant,sep=""))
    }
  }
  length(column_names)
  colnames(results_dataframe)=column_names
  
  
  
  for (tissue in unique(metadata2$tissue_type)) {
   
    #  tissue="Liver"
    #  tissue="Brain"
    
    data_trans2=as.data.frame(t(residuals_train3))
    dim(data_trans2)
    rownames(data_trans2)=metadata_train$specimen_id
    dim(data_trans2)
    
    
    # Build up-sampled model
    samplingStrategy = "up"
    rownames(metadata_train) = metadata_train$specimen_id
    TypeComparison <- metadata_train$tissue_type
    TypeString = tissue
    TypeComparisonFactor <- factor(ifelse(TypeComparison == TypeString, yes = TypeString, no = "OtherType"),
                                   levels = c(TypeString, "OtherType"))
  
    metadata_train$TypeComparison=TypeComparisonFactor
    mlDataY <- metadata_train
    mlDataX <- data_trans2[rownames(mlDataY),]
    dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
    
    trainX <- mlDataX
    trainY <- mlDataY$TypeComparison
    
    testX <- as.data.frame(t(residuals_test3))
    TypeComparison <- metadata_test$tissue_type
    TypeString = tissue
    TypeComparisonFactor <- factor(ifelse(TypeComparison == TypeString, yes = TypeString, no = "OtherType"),
                                   levels = c(TypeString, "OtherType"))
    metadata_test$TypeComparison=TypeComparisonFactor
    testY=metadata_test$TypeComparison
    
    refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
    refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
    
    refactoredTrainY <- relevel(refactoredTrainY, ref = gsub('([[:punct:]])|\\s+','',TypeString))
    refactoredTestY <- relevel(refactoredTestY, ref = gsub('([[:punct:]])|\\s+','',TypeString))
    
    #Create trainControl for Caret package
    ctrl <- trainControl(method = "repeatedcv",
                         number = 2,
                         repeats = 1,
                         summaryFunction = twoClassSummary,
                         classProbs = TRUE,
                         verboseIter = TRUE,
                         savePredictions = TRUE,
                         allowParallel=TRUE)
    
    # Build up-sampled model
    ctrl$sampling <- samplingStrategy
    print("Now training model with up sampling...")
    
    defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                                   n.trees = floor((1:3) * 50),
                                   shrinkage = 0.1,
                                   n.minobsinnode = 3)
    #Explained https://www.listendata.com/2015/07/gbm-boosted-models-tuning-parameters.html
    #Check this out https://s3.amazonaws.com/assets.datacamp.com/production/course_6650/slides/chapter2.pdf
    
    mlModel <- train(x = trainX,
                     y = refactoredTrainY,
                     method = "gbm",
                     preProcess = c("scale","center"),
                     trControl = ctrl,
                     metric = "ROC",
                     tuneGrid = defaultGBMGrid)
    
   #Calculate AUROC and AUPR
    predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[, gsub('([[:punct:]])|\\s+','',TypeString)])
    fg <- predProbs[refactoredTestY == gsub('([[:punct:]])|\\s+','',TypeString)]
    bg <- predProbs[refactoredTestY == "OtherType"]
    
    roc_GTEX=roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
    pr_GTEX=pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
    
    results_dataframe[1,paste(tissue,"_ROC",sep="")]=roc_GTEX$auc
    results_dataframe[1,paste(tissue,"_PR",sep="")]=pr_GTEX$auc.integral
    
    #Feature importance analysis
    feauter_variables= varImp(mlModel,scale=TRUE)[["importance"]]
    feauter_variables$Overall <- feauter_variables$Overall / sum(feauter_variables$Overall)
    feauter_variables$microbiomes=rownames(feauter_variables)
    feauter_variables=feauter_variables[order(feauter_variables$Overall,decreasing=TRUE),]
   
    for (contamination in c(rownames(contaminants_test2))) {
      if (contamination%in%feauter_variables$microbiomes) {
        tissue_contamination=paste(tissue,"_FIscore_",contamination,sep="")
        results_dataframe[,tissue_contamination]=feauter_variables[which(feauter_variables$microbiomes==contamination),1]
      } else {
        results_dataframe[,tissue_contamination]=0}
    }
    
  }
  
  write.csv(results_dataframe, output_file, row.names=FALSE)
}
