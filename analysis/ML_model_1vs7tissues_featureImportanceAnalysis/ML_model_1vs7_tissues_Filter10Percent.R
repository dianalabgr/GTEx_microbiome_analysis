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
if (!requireNamespace("remotes", quietly=TRUE))
  install.packages("remotes")
remotes::install_github("yiluheihei/microbiomeMarker")
library(microbiomeMarker)



#Upload metadata
metadata=read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples_QC/GTEx_Argis_QCed_all_metadata.tab")

#Upload agamenon results data

#Read agamemnon results files of all the GTEx 
GTEx_all_run = read.delim(file="../ML_model/agamemnon_after_bowtie_allSamples_species.tab")
colnames(GTEx_all_run)=gsub("\\.", "-",colnames(GTEx_all_run))
metadata$samples.submitter_id
which(colnames(GTEx_all_run)%in%metadata$samples.submitter_id)
all=GTEx_all_run[,which(colnames(GTEx_all_run)%in%c("External-ID",metadata$samples.submitter_id))]    

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

#Keep only samples from the 8 important tissues 

all_sub=all_sub[,which(tolower(metadata2$tissue_type)%in%c("heart","liver", "bladder", "muscle",
                                                           "stomach","colon","testis", "blood")),]

dim(all_sub)
dim(metadata2)
#metadata2=metadata2[which(metadata2$tissue_type%in%c("Liver","Lung","Pancreas","Brain","Heart","Bladder","Kidney","Breast","Small Intestine")),]
metadata2=metadata2[which(tolower(metadata2$tissue_type)%in%c("heart","liver", "bladder", "muscle",
                                                              "stomach","colon","testis", "blood")),]
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
for (tissue in unique(metadata$tissue_type)) {
  print(tissue)
  core_microbiome_name=paste("core_",tolower(tissue),sep="")
  core_microbiome_name=gsub(" ", "_", core_microbiome_name)
  print(core_microbiome_name)
  assign(core_microbiome_name, find_core_microbiome_Filter10(all_sub,metadata2,tissue))
  core=c(core,find_core_microbiome_Filter10(all_sub,metadata2,tissue))
}

core=unique(core)



#############################################################################################################################
################################ MACHINE LEARNING MODELS ####################################################################
#############################################################################################################################
#Parallel computing
library(parallel)
library(doMC) # for parallel computing
numCores <- detectCores()
#Set how many cores the script will use (10 cores)
registerDoMC(cores=10)


cumulative_results_ROC_PR=data.frame(matrix(nrow=length(unique(metadata2$tissue_type)), ncol=3))
colnames(cumulative_results_ROC_PR)=c("ROC","PR","PR_random")
rownames(cumulative_results_ROC_PR)=unique(metadata2$tissue_type)
for (tissue in unique(metadata2$tissue_type)) {
  #  tissue="Liver"
  #  tissue="Heart"
    
    ####################################################################################################################################
    ################################## Split the datasets to train and testing dataset 
    set.seed(42)
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
    ################################################ Separately normalise the datasets with CSS #######################################################################################
    
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
    
  
  data_trans2=as.data.frame(t(residuals_train2))
  dim(data_trans2)
  rownames(data_trans2)=metadata_train$specimen_id
  dim(data_trans2)
  
  # Build up-sampled model
  samplingStrategy = "up"
  
  rownames(metadata_train) = metadata_train$specimen_id
  TypeComparison <- metadata_train$tissue_type
  TypeString = tissue
  #  TypeString = "Liver"
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
  #  TypeString = "Liver"
  TypeComparisonFactor <- factor(ifelse(TypeComparison == TypeString, yes = TypeString, no = "OtherType"),
                                 levels = c(TypeString, "OtherType"))
  metadata_test$TypeComparison=TypeComparisonFactor
  testY=metadata_test$TypeComparison
  refactoredTrainY <- factor(gsub('([[:punct:]])|\\s+','',trainY))
  refactoredTestY <- factor(gsub('([[:punct:]])|\\s+','',testY))
  
  refactoredTrainY <- relevel(refactoredTrainY, ref = gsub('([[:punct:]])|\\s+','',TypeString))
  refactoredTestY <- relevel(refactoredTestY, ref = gsub('([[:punct:]])|\\s+','',TypeString))
  
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
  
  predProbs <- as.numeric(predict(mlModel, newdata = testX, type = "prob")[, gsub('([[:punct:]])|\\s+','',TypeString)])
  fg <- predProbs[refactoredTestY == gsub('([[:punct:]])|\\s+','',TypeString)]
  bg <- predProbs[refactoredTestY == "OtherType"]
  
  roc_GTEX=roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  pr_GTEX=pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  
  cumulative_results_ROC_PR[tissue,]=c(roc_GTEX$auc,pr_GTEX$auc.integral,pr_GTEX$rand$auc.integral)
  print(cumulative_results_ROC_PR)
  
  #Feature importance analysis
  feauter_variables= varImp(mlModel,scale=TRUE)[["importance"]]
  feauter_variables$Overall <- feauter_variables$Overall / sum(feauter_variables$Overall)
  feauter_variables$microbiomes=rownames(feauter_variables)
  feauter_variables=feauter_variables[order(feauter_variables$Overall,decreasing=TRUE),]
  write.csv(feauter_variables,file=paste("./feauture_importance_analysis/CSSOnly_FeatureImportance_Filter10Percent_",tissue,".csv",sep=""))
  
}
write.csv(cumulative_results_ROC_PR,file="./CssOnly_GTEX_results_Filter10percent_1vs7Tissues.csv")
