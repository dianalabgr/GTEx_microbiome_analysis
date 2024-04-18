`%nin%`=Negate(`%in%`)
library(parallel)
library(doMC) # for parallel computing
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

numCores <- detectCores()
#Set how many cores the script will use (10 cores)
registerDoMC(cores=10)

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

tissues=c("Heart", "Colon", "Testis", "Stomach", "Blood", "Muscle", "Liver", "Bladder")
metadata3=metadata2[which(metadata2$tissue_type%in%tissues),]
dim(metadata3)
dim(all_sub)
all_sub2=all_sub[,which(colnames(all_sub)%in%metadata3$samples.submitter_id)]
dim(all_sub2)
summary(metadata3$age_value)



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
for (tissue in tissues) {
  print(tissue)
  core_microbiome_name=paste("core_",tolower(tissue),sep="")
  core_microbiome_name=gsub(" ", "_", core_microbiome_name)
  print(core_microbiome_name)
  assign(core_microbiome_name, find_core_microbiome_Filter10(all_sub2,metadata3,tissue))
  core=c(core,find_core_microbiome_Filter10(all_sub2,metadata3,tissue))
}

core=unique(core)

####################################################################################################################################
################################## ANALYSIS PER TISSUE #############################################################################
####################################################################################################################################
cumulative_results=data.frame(matrix(nrow=1, ncol=5))
colnames(cumulative_results)=c("tissue","factor","RMSE","Rsquared","MAE")


#for (tissue in unique(metadata3$tissue_type)) {
for (tissue in tissues) {
  #tissue="Colon"
  print(tissue)
  
  #Keep the samples for the specific tissue
  dim(metadata3)
  metadata_tissue=metadata3[which(metadata3$tissue_type==tissue),]
  dim(metadata_tissue)
  
  all_sub4=all_sub2[,metadata_tissue$samples.submitter_id]
  dim(all_sub4)
  metadata4=metadata3[which(metadata3$samples.submitter_id%in%metadata_tissue$samples.submitter_id),]
  dim(metadata4)
  metadata4$samples.submitter_id==colnames(all_sub4)
  
  if (length(metadata_tissue$samples.submitter_id)<20) {
    cumulative_results=rbind(cumulative_results,data.frame("tissue"=tissue, "factor"="age", 
                                                           "RMSE"=NA, "Rsquared"=NA, "MAE"=NA))
  }else {
  ####################################################################################################################################
  ################################## Split the datasets to train and testing dataset 
  set.seed(42)
  #Create the variable to save the results 
  permutation_results=data.frame(matrix(nrow=1, ncol=4))
  colnames(permutation_results)=c("Iteration", 'RMSE', "Rsquared", "MAE")
  for (i in seq(1,100,by=1)) {
 # index <- createDataPartition(metadata4$hypertension_history, p = 0.8, list = FALSE)
  index <- createDataPartition(metadata4[,"age_value"], p = 0.7, list = FALSE)
  trainX <- all_sub4[,index]
  metadata_train <- metadata4[index,]
  testX <- all_sub4[,-index]
  metadata_test <- metadata4[-index,]
  dim(trainX)
  dim(testX)
  # table(metadata_train[,eval(disease)])
  # table(metadata_test[,eval(disease)])
  summary(metadata_train$age_value)
  summary(metadata_test$age_value)
  dim(all_sub4)
  
  #################################################################################################################################################
  ################################################ Separately normalise the datasets with CSS #######################################################################################
  
  #Train dataset
  metaSeqObject = newMRexperiment(trainX) 
  metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
  dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
  residuals_train=dge_normalised
  residuals_train[is.na(residuals_train)]=0
  colnames(residuals_train)=colnames(trainX)
  #Keep only the core microbiome
  residuals_train2=residuals_train[core,]
  #residuals_train2=residuals_train
  
  #Train dataset
  metaSeqObject = newMRexperiment(testX) 
  metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
  dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
  residuals_test=dge_normalised
  residuals_test[is.na(residuals_test)]=0
  rownames(residuals_test)
  colnames(residuals_test)=colnames(testX)
  #Keep only the core microbiome
  residuals_test2=residuals_test[core,]
  
  #residuals_test2=residuals_test
  
  
  ##################################################33333333
  ###################### Machine learning models 
  ########################################################
  #Build the model
   trainX <- as.data.frame(t(residuals_train2))
  
   testX <- as.data.frame(t(residuals_test2))
  
  dim(trainX)
  
  ctrl <- trainControl(method = "repeatedcv",
                       number = 5)
                       
  defaultGBMGrid <-  expand.grid(interaction.depth = seq(1,3),
                                 n.trees = floor((1:3) * 50),
                                 shrinkage =  0.1,
                                 n.minobsinnode = 3)
  #Explained https://www.listendata.com/2015/07/gbm-boosted-models-tuning-parameters.html
  #Check this out https://s3.amazonaws.com/assets.datacamp.com/production/course_6650/slides/chapter2.pdf
  mlModel <- train(x = trainX,
                   y = metadata_train$age_value,
                   method = "gbm",
                   trControl = ctrl,
                   tuneGrid = defaultGBMGrid)
  
  # Make predictions on the test set
  predictions <- predict(mlModel, newdata =testX)
  
  # Evaluate the model performance
  rmse <- sqrt(mean((metadata_test$age_value - predictions)^2))
  cat("Root Mean Squared Error (RMSE):", rmse, "\n")
  
  resampling_results <- data.frame(
    Observed = metadata_test$age_value,
    Predicted = predictions
  )
  
  mse <- postResample(resampling_results$Predicted, resampling_results$Observed)
  cat("Mean Squared Error (MSE):", mse, "\n")
  
  permutation_results=rbind(permutation_results, data.frame("Iteration"=i,"RMSE"=as.numeric(mse[1]),"Rsquared"=as.numeric(mse[2]),"MAE"=as.numeric(mse[3])))
  print(permutation_results)
  }
  permutation_results=permutation_results[2:101,]
  dim(permutation_results)
  write.csv(permutation_results,paste("./permutations_age/CSS_core_",sub(" ","_",tissue),"_","age",".csv",sep=""))
  temp=data.frame(tissue, "age", mean(unlist(permutation_results$RMSE)), mean(unlist(permutation_results$Rsquared)), mean(unlist(permutation_results$MAE)))
  colnames(temp) = c("tissue","factor","RMSE","Rsquared","MAE")
  cumulative_results=rbind(cumulative_results, temp)
  }
}
  cumulative_results2=cumulative_results[2:length(cumulative_results$tissue),]
  write.csv(cumulative_results2,file=paste("./CSS_total_tissues","_","age",".csv",sep=""))
  
  