#Create the database with reference and representative bacterial genomes for AGAMEMNON, following the instruction at https://github.com/ivlachos/agamemnon/wiki/Use-case


#We will use the collection of complete representative and reference genomes from NCBI RefSeq database. To retrieve these genomes in FASTA format, you need to do the #following:

#Download the assembly_summary.txt file that provides information that can be used to identify a set of genome assemblies of interest, along with their FTP file paths in the RefSeq FTP:

wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

#Extract the ftp file paths for all latest, representative and reference bacterial genomes:

awk -F "\t" '$12=="Complete Genome" && $11=="latest" && ($5 == "representative genome" || $5 == "reference genome") {print $20}' assembly_summary.txt > ftpdirpaths.txt

#Append the filename of interest (in our case _genomic.fna.gz), to each of the ftp links:

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths.txt > genomes_to_download.txt

#Create the folder "genomes" and download the genomes in there:

mkdir genomes
cd genomes
wget -i genomes_to_download.txt

#Uncompress the individual genome files and create a multi-FASTA file comprising all downloaded genomes:

gunzip *.gz
cat *.fna > microbial_reference.fa

#Optionally, to save some space, you can remove the individual .fna genome files and keep only the "microbial_reference.fa" file:

rm *.fna
