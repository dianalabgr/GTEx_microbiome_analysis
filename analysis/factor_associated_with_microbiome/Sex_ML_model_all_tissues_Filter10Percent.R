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

numCores = detectCores()
#Set how many cores the script will use (10 cores)
registerDoMC(cores=1)

#Upload metadata
metadata=read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples_QC/GTEx_Argis_QCed_all_metadata.tab")

#Upload agamenon results data

#Read agamemnon results files of all the GTEx 
GTEx_all_run = read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples/ml_models_after_bowtie/species/ML_model/agamemnon_after_bowtie_allSamples_species.tab")
colnames(GTEx_all_run)=gsub("\\.", "-",colnames(GTEx_all_run))
metadata$samples.submitter_id
which(colnames(GTEx_all_run)%in%metadata$samples.submitter_id)
all=GTEx_all_run[,which(colnames(GTEx_all_run)%in%c("External-ID",metadata$samples.submitter_id))]    


#Remove the first column of all and make the rownames of all the different species
rownames=all$"External-ID"
all=all[,2:length(colnames(all))] #Columns the samples and rows the microorganisms
rownames(all)=rownames
dim(all)

#From all remove the species with zero value in all the samples 
all_sub=all[rowSums(all==0, na.rm=TRUE)<ncol(all), ]

#Change the order of the records of metadata according to all_sub columns
metadata2=metadata[which(metadata$samples.submitter_id%in%colnames(all_sub)), ]
metadata2=metadata2[match(colnames(all_sub), metadata2$samples.submitter_id), ]
dim(metadata2)

#Keep only samples from the 8 tissues with tissue-specific microbiome
tissues=c("Heart", "Colon", "Testis", "Stomach", "Blood", "Muscle", "Liver", "Bladder")
metadata3=metadata2[which(metadata2$tissue_type%in%tissues),]
dim(metadata3)
dim(all_sub)
all_sub2=all_sub[,which(colnames(all_sub)%in%metadata3$samples.submitter_id)]
dim(all_sub2)
table(metadata3$sex)


#Find the core microbiome per tissue
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

#Find all the species which belong at least for one tissue at the core microbiome
core=unique(core)


####################################################################################################################################
################################## ANALYSIS PER TISSUE #############################################################################
####################################################################################################################################
cumulative_results=data.frame(matrix(nrow=1, ncol=7))
colnames(cumulative_results)=c("tissue","trait","ROC","margin_of_error_AUROC", "PR",
                               "margin_of_error_AUPR", "PR_random")

#for (tissue in unique(metadata3$tissue_type)) {
tissues2=tissues[which(tissues%nin%"Testis")]
for (tissue in tissues2) {
  # tissue="Heart"
  print(tissue)

  #Keep the samples for the specific tissue
  dim(metadata3)
  metadata_tissue=metadata3[which(metadata3$tissue_type==tissue),]
  dim(metadata_tissue)
  table(metadata_tissue[,"sex"])
  all_sub4=all_sub2[,metadata_tissue$samples.submitter_id]
  dim(all_sub4)
  metadata4=metadata3[which(metadata3$samples.submitter_id%in%metadata_tissue$samples.submitter_id),]
  dim(metadata4)
  metadata4$samples.submitter_id==colnames(all_sub4)
  
  if (length(metadata_tissue$samples.submitter_id)<20) {
    cumulative_results=rbind(cumulative_results,data.frame("tissue"=tissue, "trait"="sex", 
                                                           "ROC"=NA,"margin_of_error_AUROC"=NA,
                                                           "PR"=NA,"margin_of_error_AUPR"=NA, "PR_random"=NA))
  }else {
  ####################################################################################################################################
  ################################## Split the datasets to train and testing dataset 
  set.seed(42)
  #Create the variable to save the results 
  permutation_results=data.frame(matrix(nrow=1, ncol=4))
  colnames(permutation_results)=c("Iteration", 'AUROC', "AUPR", "AUPR_random")
  for (i in seq(1,100,by=1)) {
 # index = createDataPartition(metadata4$hypertension_history, p = 0.8, list = FALSE)
  index = createDataPartition(metadata4[,"sex"], p = 0.7, list = FALSE)
  trainX = all_sub4[,index]
  metadata_train = metadata4[index,]
  testX = all_sub4[,-index]
  metadata_test = metadata4[-index,]
  dim(trainX)
  dim(testX)
  table(metadata_train[,"sex"])
  table(metadata_test[,"sex"])
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
  #Keep only the species in the core 
  residuals_train2=residuals_train[core,]

  #Test dataset
  metaSeqObject = newMRexperiment(testX) 
  metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
  dge_normalised = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE))
  residuals_test=dge_normalised
  residuals_test[is.na(residuals_test)]=0
  rownames(residuals_test)
  colnames(residuals_test)=colnames(testX)
  #Keep only the species in the core 
  residuals_test2=residuals_test[core,]

  
  ##################################################33333333
  ###################### Machine learning models 
  ########################################################
  data_trans2=as.data.frame(t(residuals_train2))
  dim(data_trans2)
  # Build up-sampled model
  samplingStrategy = "up"
  
  rownames(metadata_train) = metadata_train$specimen_id
  TypeComparison = metadata_train[,eval("sex")]
  TypeString = "Male"
  TypeComparisonFactor = factor(ifelse(TypeComparison == TypeString, yes = "Male", no = "Female"),
                                 levels = c("Male", "Female"))
  table(TypeComparisonFactor)
  #  set.seed(872436)           # Set seed
  #metadata2$TypeComparison=sample(TypeComparisonFactor)
  metadata_train$TypeComparison=TypeComparisonFactor
  mlDataY = metadata_train
  mlDataX = data_trans2[rownames(mlDataY),]
  dim(mlDataY)[1] == dim(mlDataX)[1] # Sanity check
  indexSuper = 1:dim(mlDataY)[1]
  
  trainX = mlDataX
  trainY = mlDataY$TypeComparison
  
  testX = as.data.frame(t(residuals_test2))
  TypeComparison = metadata_test[,"sex"]
  TypeString = "Male"
  TypeComparisonFactor = factor(ifelse(TypeComparison == TypeString, yes = "Male", no = "Female"),
                                 levels = c("Male", "Female"))
  
  metadata_test$TypeComparison=TypeComparisonFactor
  testY=metadata_test$TypeComparison
  #length(testY)
  
  refactoredTrainY = factor(gsub('([[:punct:]])|\\s+','',trainY))
  refactoredTestY = factor(gsub('([[:punct:]])|\\s+','',testY))

  refactoredTrainY = relevel(refactoredTrainY, ref = gsub('([[:punct:]])|\\s+','',"Male"))
  refactoredTestY = relevel(refactoredTestY, ref = gsub('([[:punct:]])|\\s+','',"Male"))
  table(refactoredTestY)
  
  ctrl = trainControl(method = "repeatedcv",
                       number = 4,
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
  
  predProbs = as.numeric(predict(mlModel, newdata = testX, type = "prob")[, gsub('([[:punct:]])|\\s+','',"Male")])
  fg = predProbs[refactoredTestY == gsub('([[:punct:]])|\\s+','',"Male")]
  bg = predProbs[refactoredTestY == "Female"]
  
  roc_GTEX=roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  pr_GTEX=pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T, rand.compute=T)
  print(roc_GTEX)
  print(pr_GTEX)
  permutation_results=rbind(permutation_results, data.frame("Iteration"=i,"AUROC"=roc_GTEX$auc,"AUPR"=pr_GTEX$auc.integral,"AUPR_random"=pr_GTEX$rand$auc.integral))
  print(permutation_results)
  }
  permutation_results=permutation_results[2:101,]
  dim(permutation_results)
  write.csv(permutation_results,paste("./permutations_sex/CSS_core_",sub(" ","_",tissue),"_","sex",".csv",sep=""))
  
  # Calculate the standard error of the mean
  se_mean_X1 = sd(permutation_results$AUROC) / sqrt(length(permutation_results$AUROC))
  # Calculate the margin of error for a 95% confidence interval
  margin_of_error_X1 = qt(0.975, df = length(permutation_results$AUROC) - 1) * se_mean_X1
  
  # Calculate the standard error of the mean
  se_mean_X2 = sd(permutation_results$AUPR) / sqrt(length(permutation_results$AUPR))
  # Calculate the margin of error for a 95% confidence interval
  margin_of_error_X2 = qt(0.975, df = length(permutation_results$AUPR) - 1) * se_mean_X2
  
  temp=data.frame(tissue, "sex", mean(unlist(permutation_results$AUROC)), margin_of_error_X1,mean(unlist(permutation_results$AUPR)), margin_of_error_X2, mean(unlist(permutation_results$AUPR_random)))
  colnames(temp) = c("tissue","trait","ROC","margin_of_error_AUROC", "PR","margin_of_error_AUPR","PR_random")
  cumulative_results=rbind(cumulative_results, temp)
  }
}
cumulative_results2=cumulative_results[2:8,]
write.csv(cumulative_results2,file=paste("./CSS_total_tissues","_","sex",".csv",sep=""))

  
