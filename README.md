# GTEx_microbiome_analysis

This is the repository with all the code used in the manuscript 

### Machine learning models reveal microbial signatures in healthy human tissues, challenging the sterility of human organs  
Anargyros Skoulakis, Giorgos Skoufos, Armen Ovsepian, Artemis G. Hatzigeorgiou

This repository contains the code required to recreate the taxonomic and functional microbial profiles of samples. To achieve this, access to the GTEx data is necessary, as the accession of BAM files from RNA-seq is restricted. Follow the commands provided in the "code" folder. 
Additionally, the repository includes the code (folder "analysis") for replicating all analyses presented in the paper, utilizing the final results and metadata archived in Zenodo (link: https://zenodo.org/uploads/10980664). One can replicate the results using the final taxonomic and functional data and the available metadata of samples from Zenodo. However, the analysis aimed at distinguishing phenotypic traits and medical histories of samples from tissue microbiome requires access to the full metadata file, which is only available upon obtaining access to the GTEx data.

Specifically in the folder "code" : 
there are the scripts for creating the taxonomic and functional profiles of samples! The order of steps for recreating the taxonomic and functional profiles  is as follows: 
  1. QC_samples_subjects
  2. download_GTEx_data
  3. Quality_Trimming_reads
  4. Taxonomic_Assignment_Agamemnon (here are the folder Agamemnon_changed_folders that contains all the files different from the original GitHub repo, you just replace them in the AGAMEMNON folder - needs the extra files from the original repo to run, and the folder creation_reference_database with the instruction on how to create the database used in the study). To run AGAMEMNON see the https://github.com/ivlachos/agamemnon/wiki/Use-case
  5. Isolation_microbial_reads_remaping_bowtie
  6. humann3_after_bowtie
  7. Taxonomic_Assignment_Kaiju
  8. living_dataset (contains the pipeline for downloading and analyzing the data from the project PRJEB4337)


In the folder "analysis":
there are the scritps for replicating the results shown in the paper. For replicating the results, there is no need to recreate the taxonomic and functional profiles as they are already provided in Zenodo along with the needed GTEx samples metadata file. 
  1. number_reads_per_step (needed file from zenodo NumberReads_samples_logs.zip)
  2. alphaDiversity_and_core_microbiome
  3. ML_models_1vs27Tissues
  4. InSilicoContamination 
  5. ML_model_1vs7tissues_featureImportanceAnalysis
  6. BoxPlots_per_microbiome_8ImportantTissues
  7. ML_models_with_Kaiju
  8. ML_model_1vs7tissues_testLivingTissues
  9. Living_dataset
  10. ML_models_1vs7Tissues_GeneraLevel
  11. ML_models_1vs27Tissues_GeneraLevel
  12. ML_models_1vs7Tissues_Genes
  13. factor_associated_with_microbiome (Only for this analysis the restricted metadata file from GTEx is needed, as it is necessary to know the medical history and phenotypic trait of each subject)



