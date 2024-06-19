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
library(metagenomeSeq)


#Upload data
#As the analysis was runned in chuncks, here I upload the results of the different chunks of analysis 
#Read Kaiju results file
all=read.csv(file="./mnt/raid1/argis/GTEx/after_access/all_samples/kaiju/ML_models/CSS_only/Kaiju_taxonomic_results.csv")

#Read the phenotypes data
metadata=read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples_QC/GTEx_Argis_QCed_all_metadata.tab")
dim(metadata)

#Remove the first column of all and make the rownames of all the different species
rownames=all$taxon_name
all=all[,2:length(colnames(all))] #Columns the samples and rows the microorganisms
rownames(all)=rownames
dim(all)
rownames(all)

#Remove some non informative rows 
all=all[-which(rownames(all)%in%c("cannot be assigned to a (non-viral) species", "unclassified", "Viruses")), ]
dim(all)
#From all remove the species with zero value in all the samples 
all_sub=all[rowSums(all==0, na.rm=TRUE)<ncol(all), ]
all_sub2=all_sub 
all_sub2[is.na(all_sub2)]=0
all_sub2=all_sub2[rowSums(all_sub2==0, na.rm=TRUE)<ncol(all_sub2), ]
columns=gsub("\\.","-",colnames(all_sub2)[2:length(colnames(all_sub2))])
colnames(all_sub2)=c("taxon_name",columns)

#Change the order of the records of metadata according to all_sub columns
metadata2=metadata[match(colnames(all_sub2), metadata$specimen_id), ]
length(metadata2$specimen_id==colnames(all_sub2))


# #Keep only samples from the 8 important tissues
# 
all_sub2=all_sub2[,which(tolower(metadata2$tissue_type)%in%c("heart","liver", "bladder", "muscle",
                                                             "stomach","colon","testis", "blood")),]
dim(all_sub2)
dim(metadata2)
metadata2=metadata2[which(tolower(metadata2$tissue_type)%in%c("heart","liver", "bladder", "muscle",
                                                              "stomach","colon","testis", "blood")),]
table(metadata2$tissue_type)
dim(all_sub2)
                                                      
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
  assign(core_microbiome_name, find_core_microbiome_Filter10(all_sub2,metadata2,tissue))
  core=c(core,find_core_microbiome_Filter10(all_sub2,metadata2,tissue))
}

core=unique(core)
length(core)

#############################################################################################################################
######################################### Living DATASET ###################################################################

#Upload living dataset
data_healthy_first=read.delim("./Kaiju_concatenated_results_livingSamples.tsv")
rownames(data_healthy_first)=data_healthy_first$taxon_name
data_healthy_first=data_healthy_first[,2:length(colnames(data_healthy_first))]

#Delete the QCed_ from the names of columns of the dataset
columns_withoutQCed_healthy=unlist(lapply(colnames(data_healthy_first), function(x) substring(x,17)))
columns_withoutQCed_healthy=unlist(lapply(columns_withoutQCed_healthy, function(x) substr(x,1,nchar(x)-4)))
colnames(data_healthy_first)=columns_withoutQCed_healthy

metadata_healthy=readRDS("/mnt/raid1/argis/GTEx/after_access/normal_samples_OtherDatesets/metadata_healthy.Rdata")
dim(data_healthy_first)
which(is.na(data_healthy_first))
data_healthy_first[is.na(data_healthy_first)] = 0

#Normalise with CSS
metaSeqObject = newMRexperiment(data_healthy_first)
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject))

dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
data_normalised_healthy_first=dge_normalised
data_normalised_healthy_first[is.na(data_normalised_healthy_first)]=0
rownames(data_normalised_healthy_first)
data_normalised_healthy_first=data_normalised_healthy_first[core,]
dim(data_normalised_healthy_first)

#Check healthy metadata
metadata_healthy=metadata_healthy[which(metadata_healthy$sample%in%colnames(data_normalised_healthy_first)),]
rownames(metadata_healthy)=seq(1, length(rownames(metadata_healthy)), by=1)
table(metadata_healthy$tissue)

#Keep only the samples for the 8 tissues 
metadata_healthy=metadata_healthy[which(tolower(metadata_healthy$tissue)%in%c("heart","liver", "bladder", "muscle",
                                                                              "stomach","colon","testis", "blood")),]
data_normalised_healthy_first=data_normalised_healthy_first[,which(colnames(data_normalised_healthy_first)%in%metadata_healthy$sample)]

dim(data_normalised_healthy_first)
dim(metadata_healthy)



#############################################################################################################################
################################ MACHINE LEARNING MODELS ####################################################################
#############################################################################################################################
#Parallel computing
library(parallel)
library(doMC) # for parallel computing
numCores = detectCores()
#Set how many cores the script will use (10 cores)
registerDoMC(cores=10)

set.seed(50)
cumulative_results_ROC_PR=data.frame(matrix(nrow=length(unique(metadata2$tissue_type)), ncol=10))
colnames(cumulative_results_ROC_PR)=c("ROC","margin_of_error_ROC","PR","margin_of_error_PR","PR_random",
                                      "ROC_living_data","margin_of_error_ROC", "PR_living_data","margin_of_error_PR", "PR_living_data_random")
#rownames(cumulative_results_ROC_PR)=seq(1, 50, by=1)
rownames(cumulative_results_ROC_PR)=unique(metadata2$tissue_type)
for (tissue in unique(metadata2$tissue_type)) {
  tissue_cumulative=data.frame(matrix(nrow=100, ncol=6))
  for (i in seq(1,100,1)) {
    
    
    ####################################################################################################################################
    ################################## Split the datasets to train and testing dataset 
    index = createDataPartition(metadata2$tissue_type, p = 0.7, list = FALSE)
    trainX = all_sub2[,index]
    metadata_train = metadata2[index,]
    testX = all_sub2[,-index]
    metadata_test = metadata2[-index,]
    dim(trainX)
    dim(testX)
    
    #Train dataset
    metaSeqObject = newMRexperiment(trainX) 
    metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
    dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
    trainX_CSS=dge_normalised
    
    metaSeqObject = newMRexperiment(testX) 
    metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
    dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
    testX_CSS=dge_normalised
    
    
    ##Keep only the microbiomes in core 
    dim(t(trainX_CSS))
    residuals_train2=trainX_CSS[core,]
    dim(residuals_train2)
    residuals_test2=testX_CSS[core,]
    dim(residuals_test2)
    rownames(residuals_train2)[1:10]
    rownames(residuals_test2)[1:10]
    
  data_trans2=as.data.frame(t(residuals_train2))
  dim(data_trans2)
  # data_trans2=data_trans2[,colnames(data_trans2)%in%microbiomes]
  rownames(data_trans2)=metadata_train$specimen_id
  dim(data_trans2)
  # Build up-sampled model
  samplingStrategy = "up"
  
  rownames(metadata_train) = metadata_train$specimen_id
  TypeComparison = metadata_train$tissue_type
  TypeString = tissue
  TypeComparisonFactor = factor(ifelse(TypeComparison == TypeString, yes = TypeString, no = "OtherType"),
                                 levels = c(TypeString, "OtherType"))
  metadata_train$TypeComparison=TypeComparisonFactor
  mlDataY = metadata_train
  mlDataX = data_trans2[rownames(mlDataY),]
  dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
  
  trainX = mlDataX
  trainY = mlDataY$TypeComparison
  dim(trainX)
  
  testX = as.data.frame(t(residuals_test2))
  #dim(data_gtex_test)
  dim(testX)
  dim(metadata_test)
  TypeComparison = metadata_test$tissue_type
  TypeString = tissue
  #  TypeString = "Liver"
  TypeComparisonFactor = factor(ifelse(TypeComparison == TypeString, yes = TypeString, no = "OtherType"),
                                 levels = c(TypeString, "OtherType"))
  table(TypeComparisonFactor)
  metadata_test$TypeComparison=TypeComparisonFactor
  
  mlDataY = metadata_test
  testY = mlDataY[]$TypeComparison
  
  refactoredTrainY = factor(gsub('([[:punct:]])|\\s+','',trainY))
  refactoredTestY = factor(gsub('([[:punct:]])|\\s+','',testY))
  
  refactoredTrainY = relevel(refactoredTrainY, ref = gsub('([[:punct:]])|\\s+','',TypeString))
  refactoredTestY = relevel(refactoredTestY, ref = gsub('([[:punct:]])|\\s+','',TypeString))
  
  #Create trainControl for the Caret package
  ctrl = trainControl(method = "repeatedcv",
                       number = 2,
                       repeats = 1,
                       summaryFunction = twoClassSummary,
                       classProbs = TRUE,
                       verboseIter = TRUE,
                       savePredictions = TRUE,
                       allowParallel=TRUE)
  
  # Build up-sampled model
  ctrl$sampling = samplingStrategy
  print("Now training model with up sampling...")
  
  defaultGBMGrid =  expand.grid(interaction.depth = seq(1,3),
                                 n.trees = floor((1:3) * 50),
                                 shrinkage = 0.1,
                                 n.minobsinnode = 3)
  #Explained https://www.listendata.com/2015/07/gbm-boosted-models-tuning-parameters.html
  #Check this out https://s3.amazonaws.com/assets.datacamp.com/production/course_6650/slides/chapter2.pdf
  
  mlModel = train(x = trainX,
                   y = refactoredTrainY,
                   method = "gbm",
                   preProcess = c("scale","center"),
                   trControl = ctrl,
                   metric = "ROC",
                   tuneGrid = defaultGBMGrid)
  
  
  predProbs = as.numeric(predict(mlModel, newdata = testX, type = "prob")[, gsub('([[:punct:]])|\\s+','',TypeString)])
  fg = predProbs[refactoredTestY == gsub('([[:punct:]])|\\s+','',TypeString)]
  bg = predProbs[refactoredTestY == "OtherType"]
  
  roc_GTEX=roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  pr_GTEX=pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  print(roc_GTEX)
  print(pr_GTEX)
  
  #Check with normal living tissue from different dataset
  
  table(metadata_healthy$tissue)
  data_normalised_healthyY=metadata_healthy$tissue
  data_normalised_healthyY=factor(ifelse(tolower(data_normalised_healthyY) %in% tolower(tissue), yes = TypeString, no = "OtherType"),
                                  levels = c(TypeString, "OtherType"))
  data_normalised_healthy=as.data.frame(t(data_normalised_healthy_first))
  dim(data_normalised_healthy)
  dim(testX)
  colnames(data_normalised_healthy)
  colnames(testX)
  
  data_normalised_healthy2=data_normalised_healthy[,which(colnames(data_normalised_healthy)%in%colnames(testX))]
  data_normalised_healthy2[,colnames(testX)[which(colnames(testX)%nin%colnames(data_normalised_healthy))]]=0
  data_normalised_healthy2=data_normalised_healthy2[colnames(testX)]
  #
  dim(data_normalised_healthy2)
  dim(testX)
  colnames(data_normalised_healthy2)
  colnames(testX)
  #data_normalised_healthy=data_normalised_healthy[colnames(testX)]
  predProbs2 = as.numeric(predict(mlModel, newdata = data_normalised_healthy2, type = "prob")[, gsub('([[:punct:]])|\\s+','',TypeString)])
  length(predProbs2)
  fg2 = predProbs2[which(data_normalised_healthyY == TypeString)]
  bg2 = predProbs2[data_normalised_healthyY == "OtherType"]
  
  roc_living=roc.curve(scores.class0 = fg2, scores.class1 = bg2, curve = T)
  pr_living=pr.curve(scores.class0 = fg2, scores.class1 = bg2, rand.compute=T)
  print(roc_living)
  print(pr_living)
  tissue_cumulative[i,]=c(roc_GTEX$auc,pr_GTEX$auc.integral,pr_GTEX$rand$auc.integral,roc_living$auc,pr_living$auc.integral,pr_living$rand$auc.integral)
  print(tissue_cumulative)
  }
  write.csv(tissue_cumulative, file=paste("ROC_shuffle_",tissue,".csv",sep=""))
  
  # Calculate the standard error of the mean
  se_mean_X1 = sd(tissue_cumulative$X1) / sqrt(length(tissue_cumulative$X1))
  # Calculate the margin of error for a 95% confidence interval
  margin_of_error_X1 = qt(0.975, df = length(tissue_cumulative$X1) - 1) * se_mean_X1
  
  # Calculate the standard error of the mean
  se_mean_X2 = sd(tissue_cumulative$X2) / sqrt(length(tissue_cumulative$X2))
  # Calculate the margin of error for a 95% confidence interval
  margin_of_error_X2 = qt(0.975, df = length(tissue_cumulative$X2) - 1) * se_mean_X2
  
  # Calculate the standard error of the mean
  se_mean_X4 = sd(tissue_cumulative$X4) / sqrt(length(tissue_cumulative$X4))
  # Calculate the margin of error for a 95% confidence interval
  margin_of_error_X4 = qt(0.975, df = length(tissue_cumulative$X4) - 1) * se_mean_X4
  
  # Calculate the standard error of the mean
  se_mean_X5 = sd(tissue_cumulative$X5) / sqrt(length(tissue_cumulative$X5))
  # Calculate the margin of error for a 95% confidence interval
  margin_of_error_X5 = qt(0.975, df = length(tissue_cumulative$X5) - 1) * se_mean_X5
  
  cumulative_results_ROC_PR[tissue,]=c(mean(tissue_cumulative$X1),margin_of_error_X1,mean(tissue_cumulative$X2),margin_of_error_X2,
                                       mean(tissue_cumulative$X3),mean(tissue_cumulative$X4),margin_of_error_X4,
                                       mean(tissue_cumulative$X5),margin_of_error_X5,mean(tissue_cumulative$X6))
  print(cumulative_results_ROC_PR)

}

write.csv(cumulative_results_ROC_PR,file="./8Tissues_CssOnly_GTEX_results_Filter10percent_ALL_Kaiju_LiivingProject.csv")
