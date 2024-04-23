library(reshape2)
library(dplyr)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(data.table)

#Upload in silico contamination results
cumulative_results=read.delim(file="header_cumulative_results.csv", header=TRUE,sep=",")

#The last row contains the average of each column, I have done this on excell
mean=cumulative_results[101,]
colnames(cumulative_results)[which(colnames(cumulative_results)%like%"ROC")]
cumulative_dataframe=data.frame(matrix(nrow=12,ncol=12))
colnames(cumulative_dataframe)=c("contamination","Heart", "Liver", "Small.Intestine", "Bladder", "Brain", "Muscle", 
                                "Stomach", "Colon", "Testis", "Blood", "Salivary.Gland")                                
rownames(cumulative_dataframe)=c("contamination_1","contamination_2","contamination_3","contamination_4",
                                 "contamination_5","contamination_6","contamination_7","contamination_8",
                                 "contamination_9","contamination_10","contamination_11","contamination_12")
cumulative_dataframe$contamination=rownames(cumulative_dataframe)

#Complete the dataframe with the averages values
for (name in colnames(mean)){
  if (grepl("_FIscore_", name)) {
    tissue=unlist(strsplit(name,"_FIscore_"))[1]
    contaminant=unlist(strsplit(name,"_FIscore_"))[2]
    print(tissue)
    print(contaminant)
    cumulative_dataframe[contaminant,tissue]=mean[1,which(colnames(mean)==name)]
  }
}

# Convert dataframe to long format
cumulative_dataframe_long = melt(cumulative_dataframe, id.vars = "contamination")
cumulative_dataframe_long$value=as.numeric(cumulative_dataframe_long$value)
cumulative_dataframe_long$value = round(cumulative_dataframe_long$value*100, digits = 2)

jpeg("12_contaminants_bigger.jpg",width = 25, height = 10, units = "in", res=450)
ggplot(cumulative_dataframe_long,
      aes(x = factor(contamination,level=c("contamination_1","contamination_2","contamination_3","contamination_4",
                                           "contamination_5","contamination_6","contamination_7","contamination_8",
                                           "contamination_9","contamination_10","contamination_11","contamination_12")), y = variable, 
          size = value, label = value, fill = value)) +
  geom_point(aes(fill = value), shape = 21) + 
  geom_text(size = 7, nudge_y = -0.3) +
  scale_size(range = c(2, 15), guide = F) +
  ggpubr::rotate_x_text() +
  # coord_flip() +
  theme_pubr() + scale_fill_gradient(lim=c(0,100), low="blue",high="darkorange", name = "Percentage Contribution of Pseudo-Contaminants to Model Predictions") +
  theme(text=element_text(size=8), legend.key.size = unit(1.5, 'cm'), legend.text = element_text(size=22), legend.title = element_text(size=22)) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1, size =25), axis.text.y = element_text(size=25))  
# theme(panel.grid.major = element_line(linetype = 2, color = "black"))
dev.off()


