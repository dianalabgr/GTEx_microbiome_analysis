2 changes have been done: 
1)IN the rules folder in the script taxonomic_ranks.rule I replaced the script getTaxonomy.py with the script getTaxonomy2.py, in which in the line 93 I had added some more taxIDs in order to stop the loop. 
AND in the line 115 the funtion counts2linege has changed in order to check for the lineage of the molecules that have no zero counts.


2)IN the rules folder in the script mappings.rule I have changed the command of pufferfish. I have added in the command the argument --noOrphans in order to keep for the taxonomic assignement only the reads for which both pairs are mapped concordinately. 

3) To run Agamemnon you need to run the following command 
export LD_LIBRARY_PATH=/mnt/raid1/argis/tools/agamemnon/agamemnon/binaries/pufferfish/external/install/lib
