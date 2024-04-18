import sys
import pandas as pd
from os.path import exists

#Read kaiju result from command line
kaiju_results=pd.read_csv(sys.argv[1], sep="\t")
print(sys.argv[1])

#While if exists 
if exists("./concatenated_results.tsv") :
    concatenate=pd.read_csv("./concatenated_results.tsv", sep="\t")
    new=kaiju_results[["taxon_name","reads"]]
    name_file=kaiju_results["file"][0]
    new.rename(columns = {'reads':name_file}, inplace = True)
    concatenate_all=pd.merge(concatenate,new, how="outer", on = "taxon_name")
    concatenate_all.to_csv('./concatenated_results.tsv', sep ='\t', index=False)


#If concatenated_results.tsv does not exist, then 
if not exists("./concatenated_results.tsv") :
    print("argis")
    concatenate=kaiju_results[["taxon_name","reads"]]
    name_file=kaiju_results["file"][0]
    concatenate.rename(columns = {'reads':name_file}, inplace = True)
    concatenate.to_csv('./concatenated_results.tsv', sep ='\t', index=False)


