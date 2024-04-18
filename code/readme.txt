In this folder : 
there are the scripts for creating all the necessary files for the analysis! The order of steps for recreating the results are : 
  1. download_GTEx_data
  2. QC_samples_subjects
  3. Quality_Trimming_reads
  4. Taxonomic_Assignment_Agamemnon (here are the folder Agamemnon_changed_folders that contains all the files different from the original github repo, you just replace them in the AGAMEMNON folder - needs the extra files from the original repo to run, and the folder creation_reference_database with the instruction on how to create the database used in the study). To run AGAMEMNON see the https://github.com/ivlachos/agamemnon/wiki/Use-case
  5. Isolation_microbial_reads_remaping_bowtie
  6. humann3_after_bowtie
  7. Taxonomic_Assignment_Kaiju
  8. living_dataset (contains the pipeline for downloading and analyzing the data from the project PRJEB4337)
