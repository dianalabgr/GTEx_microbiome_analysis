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

tissues_keep=c("heart","liver","small intestine", "bladder", "brain", "muscle",
               "stomach","colon","testis", "blood","salivary gland")
all_sub=all_sub[,which(tolower(metadata2$tissue_type)%in%c("heart","liver","small intestine", "bladder", "brain", "muscle",
                                                           "stomach","colon","testis", "blood","salivary gland")),]

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
  library(metagenomeSeq)
  metaSeqObject = newMRexperiment(trainX) 
  metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
  dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
  residuals_train=dge_normalised
  residuals_train[is.na(residuals_train)]=0
  colnames(residuals_train)=colnames(trainX)
  residuals_train2=residuals_train[core,]
  
  
  #Train dataset
  library(metagenomeSeq)
  metaSeqObject = newMRexperiment(testX) 
  metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
  dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
  residuals_test=dge_normalised
  residuals_test[is.na(residuals_test)]=0
  rownames(residuals_test)
  colnames(residuals_test)=colnames(testX)
  residuals_test2=residuals_test[core,]
  
  
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
  results_dataframe=data.frame(matrix(nrow=1,ncol=length(tissues_keep)*2))
  column_names=c()
  for (tissue in unique(metadata2$tissue_type)) {
    column_names=c(column_names,paste(tissue,"_ROC",sep=""),paste(tissue,"_PR",sep=""))
    
  }
  length(column_names)
  colnames(results_dataframe)=column_names
  
  
  
  for (tissue in  unique(metadata2$tissue_type)) {
    
    #  tissue="Liver"
    #  tissue="Brain"
    
    data_trans2=as.data.frame(t(residuals_train2))
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
    
    testX <- as.data.frame(t(residuals_test2))
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
    
    #Create the trainControl for Caret package
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
    
  }
  file_output=paste("./No_contamination/No_contamination_",i,"Iteration_11Tissues.csv",sep="")
  write.csv(results_dataframe,file_output , row.names=FALSE)
}
################################################################################################################################3
######################################Compare contamination with no contamination results ######################################
###############################################################################################################################

#Create cumulative_Nocontamination dataframe
cumulative_Nocontamination=data.frame(matrix(nrow=1,ncol=22))
#Read all the no contamination results 
for (i in seq(1,100,1)) {
  print(i)
  temp=read.csv(file=paste("./No_contamination/No_contamination_",i,"Iteration_11Tissues.csv",sep=""))
  colnames(cumulative_Nocontamination)=colnames(temp)
  cumulative_Nocontamination=rbind(cumulative_Nocontamination,temp)
}
cumulative_Nocontamination=cumulative_Nocontamination[2:101,]



#Create cumulative_Contamination dataframe
cumulative_Contamination=data.frame(matrix(nrow=1,ncol=22))
#Read all the no contamination results 
for (i in seq(1,100,1)) {
  print(i)
  temp=read.csv(file=paste("../results/test_",i,"_permutations_results.log",sep=""))
  ROC_cols=which(colnames(temp)%like%c("ROC"))
  PR_cols=which(colnames(temp)%like%c("PR"))
  keep_cols=c(ROC_cols,PR_cols)   
  keep_cols=sort(keep_cols)
  temp=temp[,keep_cols]
  colnames(cumulative_Contamination)=colnames(temp)
  cumulative_Contamination=rbind(cumulative_Contamination,temp)
}
cumulative_Contamination=cumulative_Contamination[2:101,]

comparison_results=data.frame(matrix(nrow=11,ncol=8))
rownames(comparison_results)=unique(metadata2$tissue_type)
colnames(comparison_results) = c("mean_ROC_contamination", "mean_ROC_Nocontamination",
                                 "test_ROC_p.value", "test_ROC_statistic",
                                 "mean_PR_contamintaion", "mean_PR_Nocontamination",
                                 "test_PR_p.value", "test_PR_statistic")
for (tissue in rownames(comparison_results)) {
  print(tissue)
  variable_ROC=paste(gsub(" ", ".", tissue),"_ROC",sep="")
  variable_PR=paste(gsub(" ", ".", tissue),"_PR",sep="")
  test_ROC=wilcox.test(cumulative_Contamination[,variable_ROC],cumulative_Nocontamination[,variable_ROC])
  test_PR=wilcox.test(cumulative_Contamination[,variable_PR],cumulative_Nocontamination[,variable_PR])
  comparison_results[tissue,]=c(mean(cumulative_Contamination[,variable_ROC]), mean(cumulative_Nocontamination[,variable_ROC]),
                                     test_ROC$p.value, test_ROC$statistic,
                                mean(cumulative_Contamination[,variable_PR]),mean(cumulative_Nocontamination[,variable_PR]),
                                     test_PR$p.value, test_PR$statistic)
}

comparison_results
write.csv(comparison_results,"Comparisons_Contamination_vs_NoContamination.csv" , row.names=FALSE)
