#for refseq_argis3 database I have used the representatives bacterial genomes of refseq (4034) , the archaea (489) and viruses genomes of complete refseq (11259) , and the fungi (81) transcriptomes of this folder. In each folder there is a readme file for instructions how each database was created. 

cat ../Bacterial_Refseq_representatives/genomes/microbial_reference.fa ../Refseq_complete/concatenated/Archaea.fa ../Refseq_complete/concatenated/viruses.fa ./fungi/fungi.fa > Refseq_argis3.fa

#Then follow the steps of agamemnon github readme https://github.com/ivlachos/agamemnon/wiki/Use-case
