# GTEx_microbiome_analysis

This is the repository with all the code used in the manuscript 

###Machine learning models reveal microbial signatures in healthy human tissues, challenging the sterility of human organs  
Anargyros Skoulakis1,2, Giorgos Skoufos1,2, Armen Ovsepian1,2, Artemis G. Hatzigeorgiou1,2

This repo contains the code for recreating the results (to do that you should have access to the GTEx data as accession of BAM files of RNA-seq data of samples are restricted). Also, it contains the code for redoing all the analyses presented in the paper, using the final results and metadata, saved in Zenodo (link https://zenodo.org/uploads/10980664). Using the final results in Zenodo and the publicly available metadata of samples, one can replicate the results. Only the analysis for discriminating the phenotypic traits and medical history of samples demands access to the full metadata file that needs access to GTEx data. Lastly, as the publicly available metadata have different headers than the full metadata, the scripts should be modified accordingly.

In the folder code : 
there are the scripts for creating all the necessary files for the analysis! The order of steps for recreating the results is as follows: 
  1. download_GTEx_data
  2. QC_samples_subjects
  3. Quality_Trimming_reads
  4. Taxonomic_Assignment_Agamemnon (here are the folder Agamemnon_changed_folders that contains all the files different from the original GitHub repo, you just replace them in the AGAMEMNON folder - needs the extra files from the original repo to run, and the folder creation_reference_database with the instruction on how to create the database used in the study). To run AGAMEMNON see the https://github.com/ivlachos/agamemnon/wiki/Use-case
  5. Isolation_microbial_reads_remaping_bowtie
  6. humann3_after_bowtie
  7. Taxonomic_Assignment_Kaiju
  8. living_dataset (contains the pipeline for downloading and analyzing the data from the project PRJEB4337)


In the folder analysis:
there are the scritps 
  1. number_reads_per_step (needed file from zenodo NumberReads_samples_logs.zip)
  2. alphaDiversity_and_core_microbiome
  3. ML_models_1vs27Tissues
  4. InSilicoContamination
  5. ML_model_1vs7tissues_featureImportanceAnalysis
  6. BoxPlots_per_microbiome_8ImportantTissues
  7. ML_models_with_Kaiju
  8. ML_models_1vs7Tissues_GeneraLevel
  9. ML_models_1vs27Tissues_GeneraLevel
  10. ML_models_1vs7Tissues_Genes
  11. factor_associated_with_microbiome
  12. ML_model_1vs7tissues_testLivingTissues
  13. Living_dataset



