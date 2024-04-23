Here is the analysis of in silico contamination 
Run first the snakemake file snakefile_pertubations_contaminations to create the in silico contaminants 

Then create the folder results and 
Then run the R script contamination_gbm_model_contamination_afterCSS_github.R to create the results for each iteration of contamination 

The to gather all the results run the following commands in terminal

cd results

samples=($(ls))

for sample in ${samples[@]}; do echo $sample; cat $sample | tail -n +2 >> cumulative_results.csv; done

head -n 1 test_1_permutations_results.log > header.csv

cat header.csv cumulative_results.csv > header_cumulative_results.csv


