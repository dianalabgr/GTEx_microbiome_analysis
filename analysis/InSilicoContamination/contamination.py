import pandas as pd
import sys
from io import StringIO
import numpy as np
import random as random 
import sys

#Created file with filename as the first argument passed in the script
filename=sys.argv[1]

#Upload metadata 
metadata=pd.read_csv('/mnt/raid1/argis/GTEx/after_access/all_samples_QC/GTEx_Argis_QCed_all_metadata.tab',sep="\t")
print(len(metadata["date_genotype_expression"]))

#Set the seed
#random.seed(5)

#contamination_1: High-Contamination number of reads across 1000 
contamination_1_samples=random.sample(list(metadata["specimen_id"]), 1000)


#contamination_2: High-Contamination number of reads across 500 samples from 1 random BSS 
contaminated_SMCENTER=metadata["bss_collection_site"][random.randrange(len(metadata["bss_collection_site"]))]
if metadata[metadata["bss_collection_site"]==contaminated_SMCENTER].shape[0] < 500:
	contamination_2_samples=metadata[metadata["bss_collection_site"]==contaminated_SMCENTER]["specimen_id"].tolist()
else : 
	contamination_2_samples=random.sample(list(metadata[metadata["bss_collection_site"]==contaminated_SMCENTER]["specimen_id"]), 500)



#contamination_3:High-Contamination number of reads across  150 at most samples from 1 date of SMGEBTCHD (date_genotype_expression)
contaminated_SMGEBTCHD=metadata["date_genotype_expression"][random.randrange(len(metadata["date_genotype_expression"]))]
if metadata[metadata["date_genotype_expression"]==contaminated_SMGEBTCHD].shape[0] < 150:
	contamination_3_samples=metadata[metadata["date_genotype_expression"]==contaminated_SMGEBTCHD]["specimen_id"].tolist()
else : 
	contamination_3_samples=random.sample(list(metadata[metadata["date_genotype_expression"]==contaminated_SMGEBTCHD]["specimen_id"]), 150)


#Contamination_4: High-Contamination number of reads across  150 samples from 1 date of SMNABTCHD date_nucleic_acid_isolation
contaminated_SMNABTCHD=metadata["date_nucleic_acid_isolation"][random.randrange(len(metadata["date_nucleic_acid_isolation"]))]
if metadata[metadata["date_nucleic_acid_isolation"]==contaminated_SMNABTCHD].shape[0] <= 150:
	contamination_4_samples=metadata[metadata["date_nucleic_acid_isolation"]==contaminated_SMNABTCHD]["specimen_id"].tolist()
else :
	contamination_4_samples=random.sample(list(metadata[metadata["date_nucleic_acid_isolation"]==contaminated_SMNABTCHD]["specimen_id"]), 150)


#Contamination_5:High-Contamination number of reads across 150 samples at mostfrom 1 SMGEBTCH (genotype_expression_batch_id)
contaminated_SMGEBTCH=metadata["genotype_expression_batch_id"][random.randrange(len(metadata["genotype_expression_batch_id"]))]
if metadata[metadata["genotype_expression_batch_id"]==contaminated_SMGEBTCH].shape[0] <= 150:
	contamination_5_samples=metadata[metadata["genotype_expression_batch_id"]==contaminated_SMGEBTCH]["specimen_id"].tolist()
else :
	contamination_5_samples=random.sample(list(metadata[metadata["genotype_expression_batch_id"]==contaminated_SMGEBTCH]["specimen_id"]), 150)


#Contamintaion_6:High-Contamination number of reads across 150 samples at most, from 1 random SMNABTCH (nucleic_acid_isolation_batch_id)
contaminated_SMNABTCH=metadata["nucleic_acid_isolation_batch_id"][random.randrange(len(metadata["nucleic_acid_isolation_batch_id"]))]
if metadata[metadata["nucleic_acid_isolation_batch_id"]==contaminated_SMNABTCH].shape[0] <= 150:
	contamination_6_samples=metadata[metadata["nucleic_acid_isolation_batch_id"]==contaminated_SMNABTCH]["specimen_id"].tolist()
else :
	contamination_6_samples=random.sample(list(metadata[metadata["nucleic_acid_isolation_batch_id"]==contaminated_SMNABTCH]["specimen_id"]), 150)


#Contamination_7: Low-Contamination number of reads across all samples 
contamination_7_samples=metadata["specimen_id"].tolist()
print("contamination 7"+str(len(contamination_7_samples)))


# Contamination_8: Low-Contamination number of reads across all samples from 1 random BSS
contaminated_8_SMCENTER=metadata["bss_collection_site"][random.randrange(len(metadata["bss_collection_site"]))]
contamination_8_samples=metadata[metadata["bss_collection_site"]==contaminated_8_SMCENTER]["specimen_id"].tolist()


#Contamination_9: Low-Contamination number of reads across all samples from 1 date of SMGEBTCHD (date_genotype_expression)
contaminated_9_SMGEBTCHD=metadata["date_genotype_expression"][random.randrange(len(metadata["date_genotype_expression"]))]
contamination_9_samples=metadata[metadata["date_genotype_expression"]==contaminated_9_SMGEBTCHD]["specimen_id"].tolist()


#Contamination_10 : Low-Contamination number of reads across all samples from 1 date of SMNABTCHD (date_nucleic_acid_isolation)
contaminated_10_SMNABTCHD=metadata["date_nucleic_acid_isolation"][random.randrange(len(metadata["date_nucleic_acid_isolation"]))]
contamination_10_samples=metadata[metadata["date_nucleic_acid_isolation"]==contaminated_10_SMNABTCHD]["specimen_id"].tolist()


# Contamination_11 : Low-Contamination number of reads across all samples from 1 random batch of SMGEBTCH (genotype_expression_batch_id)
contaminated_11_SMGEBTCH=metadata["genotype_expression_batch_id"][random.randrange(len(metadata["genotype_expression_batch_id"]))]
contamination_11_samples=metadata[metadata["genotype_expression_batch_id"]==contaminated_11_SMGEBTCH]["specimen_id"].tolist()


#Contamination_12 : Low-Contamination number of reads across all samples from 1 random batch of SMNABTCH (nucleic_acid_isolation_batch_id)
contaminated_12_SMNABTCH=metadata["nucleic_acid_isolation_batch_id"][random.randrange(len(metadata["nucleic_acid_isolation_batch_id"]))]
contamination_12_samples=metadata[metadata["nucleic_acid_isolation_batch_id"]==contaminated_12_SMNABTCH]["specimen_id"].tolist()


contaminated= [contamination_1_samples, contamination_2_samples, contamination_3_samples, contamination_4_samples, 
				contamination_5_samples, contamination_6_samples, contamination_7_samples, contamination_8_samples, 
				contamination_9_samples, contamination_10_samples, contamination_11_samples, contamination_12_samples]


contaminated_dataFrame=pd.DataFrame({'Contamination':[sublist for sublist in contaminated]})

contaminated_dataFrame.to_csv(filename, index=True)



