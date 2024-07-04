In this folder : 
there are the scripts for QC of samples, downloading GTEx data, QC of the reads, and taxonomic and functional profiling of the samples! The order of steps for recreating the results are : 
  1. QC_samples_subjects
  2. download_GTEx_data
  3. Quality_Trimming_reads
  4. Taxonomic_Assignment_Agamemnon (here are the folder Agamemnon_changed_folders that contains all the files different from the original github repo (https://github.com/ivlachos/agamemnon), you just replace them in the AGAMEMNON folder and you follow the instructions in the folder creation_reference_database on how to create the database used in the study. To run AGAMEMNON see the https://github.com/ivlachos/agamemnon/wiki/Use-case
  5. Isolation_microbial_reads_remaping_bowtie
  6. humann3_after_bowtie
  7. Taxonomic_Assignment_Kaiju
  8. living_dataset (contains the pipeline for downloading and analyzing the data from the NCBI Bioproject PRJEB4337)
