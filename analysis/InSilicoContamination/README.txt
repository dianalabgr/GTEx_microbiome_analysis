cd results

samples=($(ls))

for sample in ${samples[@]}; do echo $sample; cat $sample | tail -n +2 >> cumulative_results.csv; done

head -n 1 test_1_permutations_results.log > header.csv

cat header.csv cumulative_results.csv > header_cumulative_results.csv


