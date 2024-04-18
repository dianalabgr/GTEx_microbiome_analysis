#Create the database for reference using all the complete viruses genomes from Refseq. insight https://www.nature.com/articles/s41596-022-00738-y

wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt

awk -F "\t" '$12=="Complete Genome" && $11=="latest" {print $20}' assembly_summary.txt > ftpdirpaths.txt

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths.txt > genomes_to_download.txt

mkdir genomes
cd genomes
wget -i ../genomes_to_download.txt
