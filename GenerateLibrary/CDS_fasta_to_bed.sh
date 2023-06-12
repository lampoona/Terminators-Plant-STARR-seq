cd ~/Terminators/Trial/FinalTest/

module load seqtk/1.3

seqtk seq -F '#' TAIR_CDS_589_sample.fa > TAIR_CDS_589_sample.fastq

seqtk seq -F '#' Maize_CDS_589_sample.fa > Maize_CDS_589_sample.fastq

bowtie2 -p 5 --xeq --no-mixed -x /net/gs/vol1/home/kbubb/queitschlab/TAIR_RELEASE_GENOMES/TAIR10/TAIR10.fas -U TAIR_CDS_589_sample.fastq 2> alignment_TAIR_CDS_589_sample.txt | samtools view -bS - > TAIR_CDS_589_sample.bam  

bwa mem -t 4 -R "@RG\tID:A8100\tSM:A8100" /net/gs/vol1/home/kbubb/queitschlab2/Maize/Zea_mays.AGPv4.chr0.fa Maize_CDS_589_sample.fastq 2> alignment_Maize_CDS_589_sample.txt | samtools view -bS - > Maize_CDS_589_sample.bam
 
bedtools bamtobed -i TAIR_CDS_589_sample.bam > TAIR_CDS_589_sample.bed

bedtools bamtobed -i Maize_CDS_589_sample.bam > Maize_CDS_589_sample.bed
