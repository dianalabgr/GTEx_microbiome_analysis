Configure the profile for gen3-client
../../tools/gen3_client/gen3-client configure --profile=GTEX_microbiome --cred=./credentials_2.json --apiendpoint=https://gen3.theanvil.io

The file files_objectID_all.txt was taken from the file sequencing.tsv and contains all the objectIDs from the column object ID and the sample names from the column samples.submitter.id. (NEEDED in the pipeline). The 2 columsn are separated by space and there is no extra lines or headers.

