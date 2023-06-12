############### Arabidopsis ######################
module load bedtools/2.25.0 bedops/2.4.35 python/3.7.7 bedops/2.4.35 pandas/1.3.1
module load pyfaidx/0.5.9.2
module load pcre2/10.35 hdf5/1.10.1 R/4.0.0


wget "https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff" 
grep "three_prime_UTR" TAIR10_GFF3_genes.gff | awk -v OFS='\t' '{print $1,$4,$5,$9,$8,$7}' | awk -vOFS="\t" '{sub("Parent=","",$4); print}' | grep '\.1' --no-group-separator | awk '{if($2 != $3) print $0}' | sort-bed - > TAIR10_GFF3_genes_three_prime_UTR.gff 

#get rid of duplicates by name (sort bed does it by start/chromosome)
Rscript Dedup.R 

python last_exon.py TAIR10_GFF3_genes_three_prime_UTR.bed

#go 150 back and 20 forward
awk -v OFS='\t' '{if($6 == "+") print $1,$3-150,$3+20,$4,$5,$6; if($6 =="-") print $1,$2-20,$2+150,$4,$5,$6}' TAIR10_3UTR_last_exon.bed > Arabidopsis_170.bed 

cp /net/gs/vol1/home/sgorji/Terminators/Thomas2012_SuppTables.xlsx . 

Rscript Arabidopsis.R

bedtools getfasta -fi /net/gs/vol1/home/kbubb/queitschlab/TAIR_RELEASE_GENOMES/TAIR10/TAIR10.fas.fa -bed Arabidopsis_PACs_Thomas_plus_TAIR.bed  -s -name -fo Arabidopsis_PACs_Thomas_plus_TAIR.fa

awk '{printf "%s%s",$0,NR%2?"\t":RS}' Arabidopsis_PACs_Thomas_plus_TAIR.fa > R_Arabidopsis_PACs_Thomas_plus_TAIR.fa

Rscript Dedup_seq.R


bedtools getfasta -fi /net/gs/vol1/home/kbubb/queitschlab/TAIR_RELEASE_GENOMES/TAIR10/TAIR10.fas.fa -bed Arabidopsis_PACs_Thomas_plus_TAIR.bed  -s -name -fo Arabidopsis_PACs_Thomas_plus_TAIR.fa

python3 oligo_frequency.py Arabidopsis_PACs_Thomas_plus_TAIR.fa 170

######################## MAIZE #####################

cp /net/gs/vol1/home/sgorji/Terminators/File_S3_Maize_clean.txt . 

Rscript Maize.R 

bedtools getfasta -fi /net/gs/vol1/home/kbubb/queitschlab2/Maize/Zea_mays.AGPv4.chr0.fa -bed Maize_Top_PACs_Final.bed -s -name -fo Maize_Top_PACs_Final.fa

python3 oligo_frequency.py Maize_Top_PACs_Final.fa 170



########get RE's removed 

./Remove_RE_add_adapter.sh 

#############################################
#CONTROLS
############################################# 



#######CDS controls Arabidopsis 

Rscript Clean_TAIR_CDS.R

grep '\.1' --no-group-separator TAIR10_CDS_170_plus.bed > TAIR10_CDS_170_plus_1.bed

bedtools getfasta -fi /net/gs/vol1/home/kbubb/queitschlab/TAIR_RELEASE_GENOMES/TAIR10/TAIR10.fas.fa -bed  TAIR10_CDS_170_plus_1.bed  -s -name -fo TAIR10_CDS_170_plus.fa

awk '{printf "%s%s",$0,NR%2?"\t":RS}' TAIR10_CDS_170_plus.fa > RR_TAIR10_CDS_170_plus.fa

Rscript TAIR_CDS_Sample.R


########CDS control maize 


#V4 
curl http://ftp.ebi.ac.uk/ensemblgenomes/pub/release-39/plants/gff3/zea_mays/Zea_mays.AGPv4.39.chr.gff3.gz | unpigz > Maize_annotation_v4.gff3
grep "CDS" Maize_annotation_v4.gff3 | grep "_T001" | awk '$1="chr0"$1' | awk -v OFS="\t" '{print $1,$4,$5,$9,$6,$7}' | awk -vOFS="\t" '{split($4,a,";"); print $1,$2,$3,a[2],$5,$6}' | awk -vOFS="\t" '{sub("Parent=transcript:", "", $4); print}' | awk '{if($2 != $3) print $0}' | sort-bed - > maize_V4_CDS_clean.bed
sed -i 's/chr010/chr10/g' maize_V4_CDS_clean.bed
sed -i '/chr0Mt/d;/chr0Pt/d' maize_V4_CDS_clean.bed

sort-bed maize_V4_CDS_clean.bed > maize_V4_CDS_clean_sorted.bed

Rscript Clean_Maize_CDS.R

bedtools getfasta -fi /net/gs/vol1/home/kbubb/queitschlab2/Maize/Zea_mays.AGPv4.chr0.fa -bed maize_V4_CDS_170_length.bed -s -name -fo maize_V4_CDS_170_length.fa
awk '{printf "%s%s",$0,NR%2?"\t":RS}' maize_V4_CDS_170_length.fa > RR_maize_V4_CDS_170_length.fa

Rscript Maize_CDS_Sample.R


################################################
#randomly generated contorls 
################################################

#python script generates random oligos of length 170 with desired GC content  
#edited to take out RE sites 

python3 Generate_oligo_by_GC.py .3 170 200 > GC_30_170.fa
python3 Generate_oligo_by_GC.py .4 170 200 > GC_40_170.fa
python3 Generate_oligo_by_GC.py .5 170 200 > GC_50_170.fa
python3 Generate_oligo_by_GC.py .6 170 200 > GC_60_170.fa
python3 Generate_oligo_by_GC.py .7 170 200 > GC_70_170.fa


#python script that generates a position specific nucleotide composition sequence
python3 Position_specific_GC.py Arabidopsis_PACs_Thomas_plus_TAIR_oligo_frequency_table.txt 1000 > Arabidopsis_PACs_Thomas_plus_TAIR_Positional_random.fa
python3 Position_specific_GC.py Maize_Top_PACs_Final_oligo_frequency_table.txt 1000 > Maize_Top_PACs_Final_Positional_random.fa

python3 oligo_frequency.py Arabidopsis_PACs_Thomas_plus_TAIR_Positional_random.fa 170
python3 oligo_frequency.py Maize_Top_PACs_Final_Positional_random.fa 170


#Generate random control sequences that have the same overall ACGT percentage as the FASTA file given to the file#' 
#this version takes out any sequences with RE sites 
python3 Generate_oligo_by_seq.py Maize_Top_PACs_Final.fa 170 200 Maize > Maize_ACGT_random.fa 
python3 Generate_oligo_by_seq.py Arabidopsis_PACs_Thomas_plus_TAIR.fa 170 200 Arabidopsis > Arabidopsis_ACGT_random.fa

################################################
#G4 additional analysis 
################################################

python g4scan.py Maize_Top_PACs_Final.fa > Maize_Top_PACs_G4.bed
python g4scan.py Arabidopsis_PACs_Thomas_plus_TAIR.fa > Arabidopsis_PACs_Thomas_plus_TAIR_G4.bed


awk '{printf "%s%s",$0,NR%2?"\t":RS}' Maize_Top_PACs_Final_Final.fa > RR_Maize_Top_PACs_Final_Final.fa
awk '{printf "%s%s",$0,NR%2?"\t":RS}' Arabidopsis_PACs_Thomas_plus_TAIR_final.fa > RR_Arabidopsis_PACs_Thomas_plus_TAIR_final.fa

Rscript G4.R 


cat Maize_AAAA_G4.fa Maize_BBBB_G4.fa Maize_ABBB_G4.fa Maize_AABB_G4.fa > Maize_G4_combined.fa 
cat Arab_AAAA_G4.fa Arab_BBBB_G4.fa Arab_ABBB_G4.fa Arab_AABB_G4.fa > Arab_G4_combined.fa


#pool all controls together 
cat TAIR_CDS_589_sample.fa Maize_CDS_589_sample.fa  GC_30_170.fa GC_40_170.fa GC_50_170.fa GC_60_170.fa GC_70_170.fa Arabidopsis_PACs_Thomas_plus_TAIR_Positional_random.fa Maize_Top_PACs_Final_Positional_random.fa Maize_ACGT_random.fa Arabidopsis_ACGT_random.fa Maize_G4_combined.fa Arab_G4_combined.fa > controls.fa

awk -v OFS='\t' '{
    NAME=$1
    getline
    SEQ="TAAG"$1"AGGT"
    MUT=""
    if(SEQ ~ /((G|C)GTCTC)|(GAGAC(C|G))/) {
      while(index(SEQ, "GGTCTC") > 0) {MUT=MUT "T" index(SEQ, "GGTCTC") -2 ">A;"; sub(/GGTCTC/, "GGACTC", SEQ)};
      while(index(SEQ, "CGTCTC") > 0) {MUT=MUT "T" index(SEQ, "CGTCTC") - 2 ">A;"; sub(/CGTCTC/, "CGACTC", SEQ)};
      while(index(SEQ, "GAGACC") > 0) {MUT=MUT "A" index(SEQ, "GAGACC") - 1 ">T;"; sub(/GAGACC/, "GAGTCC", SEQ)};
      while(index(SEQ, "GAGACG") > 0) {MUT=MUT "A" index(SEQ, "GAGACG") - 1 ">T;"; sub(/GAGACG/, "GAGTCC", SEQ)};
      print substr(NAME, 2), substr(MUT, 1, length(MUT) - 1) > "RE_mutations_controls.tsv"
    }
    SEQ="GCGCCGTCTCC"SEQ"CGAGACGGTGC"
    print NAME"\n"SEQ
  }' controls.fa \
  > controls_Final.fa

#cat Maize_Top_PACs_Final_Final.fa Arabidopsis_PACs_Thomas_plus_TAIR_final.fa > Terminators.fa 

cat Maize_Top_PACs_Final_Final.fa Arabidopsis_PACs_Thomas_plus_TAIR_final.fa | awk -v OFS="\t" '{if(NR %2 ==1) name=$1; else if ($1 ~ /^[ACGT]+$/) print name"\n"$1}' > Terminators.fa 

cat Terminators.fa controls_Final.fa > Final_tmp.fa

awk '{printf "%s%s",$0,NR%2?"\t":RS}' Final_tmp.fa > RR_Final_terminator_sequences.fa

Rscript Dedup_final.R



