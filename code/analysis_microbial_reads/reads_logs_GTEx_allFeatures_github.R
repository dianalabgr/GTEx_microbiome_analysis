library(ggplot2)
library(stringr)
library(viridis)

#Upload metadata
metadata=read.delim(file="/mnt/raid1/argis/GTEx/after_access/all_samples_QC/GTEx_Argis_QCed_all_metadata.tab")

#Finding the number of total reads per sample from the metadata file
total_reads=metadata[,c("tissue_type","tissue_type_detail","specimen_id","total_reads","failed_vendor_qc_check")]
total_reads$BAM_reads=total_reads$total_reads+total_reads$failed_vendor_qc_check

#From the demoRun I used the file DemoRun_total_unmapped.log, but DemoRun contains also
#some samples than not passed the QC, so I will remove them
DemoRun_total_unmapped=read.table("DemoRun_total_unmapped.log")
DemoRun_total_unmapped=DemoRun_total_unmapped[DemoRun_total_unmapped$V1%in%metadata$samples.submitter_id,]
#From the folder logs_unmapped the file 11822samples_total_unmapped_log was used
samples11822_total_unmapped=read.table("11822samples_total_unmapped.log")
#Find if there are any common samples
which(samples11822_total_unmapped$V1%in%DemoRun_total_unmapped$V1)
samples11822_total_unmapped[2362,]
DemoRun_total_unmapped[which(DemoRun_total_unmapped$V1%in%samples11822_total_unmapped$V1),]
#there is one, so I remove it and I concatenate the files
total_unmapped=rbind(samples11822_total_unmapped[-2362,],DemoRun_total_unmapped)
colnames(total_unmapped)= c("SAMPID","unmapped_reads")
  
#Find the QCed reads per sample
#From the demoRun I used the file DemoRun_total_QCed.log, but DemoRun contains also
#some samples than not passed the QC, so I will remove them
DemoRun_total_QCed=read.table("DemoRun_total_QCed.log")
DemoRun_total_QCed=DemoRun_total_QCed[DemoRun_total_QCed$V1%in%metadata$samples.submitter_id,]
#From the folder logs_QCed the file 11822samples_total_QCed_log was used
samples11822_total_QCed=read.table("11822samples_total_QCed.log")
#Find if there are any common samples
which(samples11822_total_QCed$V1%in%DemoRun_total_QCed$V1)
samples11822_total_QCed[2362,]
DemoRun_total_QCed[which(DemoRun_total_QCed$V1%in%samples11822_total_QCed$V1),]
#there is one, so I remove it and I concatenate the files
total_QCed=rbind(samples11822_total_QCed[-2362,],DemoRun_total_QCed)
colnames(total_QCed)= c("SAMPID","QCed_reads")

#Reads after finding only the reads that are used by pufferfish
DemoRun_total_humann=read.table("super_demo_run_total_microbial_logs.log")
DemoRun_total_humann=DemoRun_total_humann[DemoRun_total_humann$V1%in%metadata$samples.submitter_id,]
samples11822_total_humann=read.table("all_samples_total_microbial_logs.log")
#Find if there are any common samples
which(samples11822_total_humann$V1%in%DemoRun_total_humann$V1)
samples11822_total_humann[2281,]
DemoRun_total_humann[which(DemoRun_total_humann$V1%in%samples11822_total_humann$V1),]
#there is one, so I remove it and I concatenate the files
total_humann=rbind(samples11822_total_humann[-2281,],DemoRun_total_humann)
colnames(total_humann)= c("SAMPID","humann_reads")

#Reads after bowtie 
bowtie_reads=read.table("after_bowtie_all_samples.log")
bowtie_reads_correct=data.frame(V1=bowtie_reads$V2,V2= bowtie_reads$V1)
bowtie_reads_correct$V1=str_sub(bowtie_reads_correct$V1,1,-12)
total_bowtie=bowtie_reads_correct[bowtie_reads_correct$V1%in%metadata$samples.submitter_id,]
colnames(total_bowtie)= c("SAMPID","unmapped_bowtie_reads")

#Microbial Reads after Agamemnon
agamemnon=read.delim("agamemnon_after_bowtie_allSamples_species.tab")
colnames(agamemnon)=gsub("\\.", "-",colnames(agamemnon))
agamemnon2=agamemnon[,colnames(agamemnon)%in%metadata$samples.submitter_id]
total_microbial=data.frame(V1=colnames(agamemnon2), V2=colSums(agamemnon2)*2)
colnames(total_microbial)= c("SAMPID","microbial_reads")


total_merge1=merge(total_reads,total_unmapped,all=TRUE, by.x="specimen_id" , by.y="SAMPID")
total_merge2=merge(total_merge1,total_QCed,all=TRUE, by.x="specimen_id" , by.y="SAMPID")
total_merge2_2=merge(total_merge2,total_humann,all=TRUE, by.x="specimen_id" , by.y="SAMPID")
total_merge3=merge(total_merge2_2,total_bowtie,all=TRUE, by.x="specimen_id" , by.y="SAMPID")
total_merge4=merge(total_merge3,total_microbial,all=TRUE, by.x="specimen_id" , by.y="SAMPID")

total_means=data.frame(matrix(ncol=9))
#total_means=data.frame(matrix(ncol=7))


#WWith tissue general
total_means=data.frame(matrix(ncol=13))
#total_means=data.frame(matrix(ncol=7))
for (tissue in unique(total_merge4$tissue_type)) {
  print(tissue)
  BAM_means=mean(total_merge4[which(total_merge4$tissue_type==tissue),"BAM_reads"])
  log_BAM_means=log10(BAM_means)
  Unmapped_means=mean(total_merge4[which(total_merge4$tissue_type==tissue),"unmapped_reads"])
  log_Unmapped_means=log10(Unmapped_means)
  QCed_means=mean(total_merge3[which(total_merge3$tissue_type==tissue),"QCed_reads"])
  log_QCed_means=log10(QCed_means)
  Microbial_means=mean(total_merge4[which(total_merge4$tissue_type==tissue),"humann_reads"])
  log_Microbial_means=log10(Microbial_means)
  after_bowtie_means=mean(total_merge4[which(total_merge4$tissue_type==tissue),"unmapped_bowtie_reads"])
  log_after_bowtie_means=log10(after_bowtie_means)
  classified_microbial_means=mean(total_merge4[which(total_merge4$tissue_type==tissue),"microbial_reads"])
  log_classified_microbial_means=log10(classified_microbial_means)
  #total_means[nrow(total_means)+1,]=c(tissue,BAM_means,log_BAM_means,Unmapped_means,log_Unmapped_means,
  #                                    Microbial_means,log_Microbial_means)
  total_means[nrow(total_means)+1,]=c(tissue,BAM_means,log_BAM_means,Unmapped_means,log_Unmapped_means,
                                      QCed_means,log_QCed_means,Microbial_means,log_Microbial_means,after_bowtie_means,log_after_bowtie_means,
                                      classified_microbial_means,log_classified_microbial_means)
}

#total_means=total_means[2:9,]  
colnames(total_means)=c("tissue", "Total_Reads", "Log Total_Reads","Unmapped_reads","Log Unmapped_reads",
                        "QCed_reads",  "Log QCed_reads", "Microbial_reads", "Log Microbial_reads", 
                        "After Bowie Microbial reads", "Log Non Human Microbial reads", "Classified Microbial Reads",
                        "Log Classified Microbial Reads")
total_means=total_means[2:29,]  
#colnames(total_means)=c("tissue", "Total_Reads", "Log Total_Reads","Unmapped_reads","Log Unmapped_reads",
#                        "Microbial_reads",  "Log Microbial_reads")

#Found the average per sample values 
#Total_Reads
summary(as.numeric(total_means[,2]))
#Unmapped Reads
summary(as.numeric(total_means[,4]))
#QCed reads
summary(as.numeric(total_means[,6]))
#Microbial Reads
summary(as.numeric(total_means[,8]))
#After bowtie non human reads
summary(as.numeric(total_means[,10]))
#Classified reads 
summary(as.numeric(total_means[,12]))


library(reshape2)
dfm= melt(total_means[,c("tissue","Log Total_Reads","Log Unmapped_reads","Log QCed_reads","Log Microbial_reads","Log Non Human Microbial reads","Log Classified Microbial Reads")], id.vars=1)
#dfm= melt(total_means[,c("tissue","Log Total_Reads","Log Unmapped_reads","Log Microbial_reads")], id.vars=1)
dfm$value=as.numeric(dfm$value)
ggplot(dfm,aes(x = tissue,y = value)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge") + 
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,10)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
jpeg(file="Reads_per_tissue_log_all_github.jpeg", width=20, height=10, units="in", res=600)
ggplot(dfm,aes(x = tissue,y = value)) + 
  geom_bar(aes(fill = variable),stat = "identity",position = "dodge") + 
  scale_fill_brewer(palette = "Dark2")+
  labs(x = "", y = "log(Number of Reads)")+
  scale_y_continuous(breaks=c(0,2,4,6,8,10), limits=c(0,9)) + theme_bw()+ 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  labs(y="log10(Number of Reads)") + theme(text = element_text(size = 28))
dev.off()





