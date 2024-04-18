Created at 15-02-2023

#Downlaod the assembly summary of Refseq for fungi
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt

#keep only the latest assemblies (but keep complete, chromosome, scaffold and contig level assemblies)
awk -F "\t" '($12=="Complete Genome" || $12 =="Chromosome") && $11=="latest" && ($5 == "representative genome" || $5 == "reference genome") {print $20}' assembly_summary.txt  > ftpdirpaths.txt
#81 transcriptomes will be downloaded 

#Create the adreeses in order to download the transcriptome for each assembly. Take only the rna.fna.gz record for each assembly
 awk 'BEGIN{FS=OFS="/";filesuffix="rna.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths.txt > genomes_to_download.txt

#Download all the rna.fasta for all assemblies.
mkdir genomes 
cd genomes
wget -i ../genomes_to_downlaod.txt

gunzip *.gz
cat *.fna > ../fungi.fa
