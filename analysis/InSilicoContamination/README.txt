Here is the analysis of in silico contamination 
1. Run first the snakemake file snakefile_pertubations_contaminations to create the in silico contaminants 
snakemake --snakefile snakefile_pertubations_contaminations --cores 4

2. Then create the folder results and 
3. Then run the R script contamination_gbm_model_contamination_afterCSS_github.R to create the results for each iteration of contamination 

4. Then to gather all the results run the following commands in terminal

cd results

samples=($(ls))

for sample in ${samples[@]}; do echo $sample; cat $sample | tail -n +2 >> cumulative_results.csv; done

head -n 1 test_1_permutations_results.log > header.csv

cat header.csv cumulative_results.csv > header_cumulative_results.csv


5. Run the script noContamination/No_contamination_gbm_model_contamination_afterCSS_github.R for creating 100 iteration without contaminations in order to compare the performances 

6. Run the R script buble_plots_github.R in order to reproduce the bubble plots 
