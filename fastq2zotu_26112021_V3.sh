#>>>>>>>>>>By Krych & Castro-Mejía 2021 <<<<<<<<<<<<<<
#>>>>>>>>>>License : University of Copenhagen<
#>>>>>>>>>>krych@food.ku.dk, jcame@food.ku.dk<<<<<<<<<<<<<<<<<<<


    #Adjusting file labels, relabeling sequences & merging_ovelaping_pairs##############################################################################

reverse_file=(*R2*.fastq)
for file1 in *R1*.fastq
do

string=${file1}
string2=$(sed "-es/_R[0-9].*\.fastq//; s/combined_//" <<< $string)
string3=$(sed "-es/-/_/g" <<< $string2)
echo "sampleID: ${string3}"

string4=$(sed "-es/S_[0-9].//g" <<< $string3)
echo ""
echo "project name: ${string4}"
echo ""

usearch -fastq_mergepairs ${file1} -reverse ${reverse_file[i++]} -fastqout merged-${file1}.fq -fastq_minovlen 100 -relabel ${string4}.

done

    #Cleaning - moving fastq files to a storage directory##############################################################################

mkdir ${string4}
mv *.fastq ${string4}/

    #Quality check of merged sequences##############################################################################

for file in merged*fq
do
usearch -fastq_filter $file -fastq_maxee 2.0  -fastq_truncqual 4 -fastq_minlen 130 -fastaout ${file}.fna
done

    #Concatenating high quality merged reads

cat merged*.fna > qual-filt.fna
rm merged*.fna
rm merged*.fq

    #Finding unique sequences##############################################################################

usearch -fastx_uniques qual-filt.fna -sizeout -fastaout uniques.fa

    #Identifying zero-radious OTUs and removing chimeric sequences

usearch -unoise3 uniques.fa -zotus zotus.fa
sed "s/Zotu/zOTU_/g" zotus.fa > zotus2.fa 
rm zotus.fa
rm uniques.fa

    #Building an OTU table##############################################################################

usearch -otutab qual-filt.fna -otus zotus2.fa -otutabout raw_zOTU.txt  -id 0.999
cat raw_zOTU.txt | sed 's/#OTU ID/OTUID/g' > raw_zOTU2.txt 
rm raw_zOTU.txt
rm qual-filt.fna

    #Taxonomic affiliation of zOTUs##############################################################################

usearch -sintax zotus2.fa -db Databases_16s/Ez-taxon-format-sintax.fasta  -strand both -tabbedout sintaxEt.txt -sintax_cutoff 0.8
cat sintaxEt.txt | cut -f1,4 | sed 's/d:/k__/g' | sed 's/p:/p__/g' | sed 's/c:/c__/g' | sed 's/o:/o__/g' | sed 's/f:/f__/g' | sed 's/g:/g__/g' | sed 's/s:/s__/g' | sed 's/,/;/g' > Etsintax.txt 

rm sintaxEt.txt
TAB=$'\t'

echo 'OTUID'"${TAB}"'taxonomy' > headers.txt 
cat headers.txt Etsintax.txt  > EtsintaxR.txt ; rm  headers.txt ; rm Etsintax.txt 


echo '#!/usr/bin/env Rscript
X = read.table("raw_zOTU2.txt", header = TRUE, sep = "\t")
Y <- read.table("EtsintaxR.txt", header = TRUE, sep = "\t", na.string = "")
table <- merge(X,Y, by = "OTUID")
write.table(table, "OTU_completed_EZtaxon.txt", sep="\t", row.names = F, quote =F)' > biom.R

chmod 755 biom.R
./biom.R
rm biom.R
rm raw_zOTU2.txt

mkdir Results_${string4}
mkdir Results_${string4}/OTU-tables
mv OTU_completed_EZtaxon.txt Results_${string4}/OTU-tables/

    #Generating phylogenetic tree ##############################################################################

mkdir Results_${string4}/trees
usearch -calc_distmx zotus2.fa -tabbedout mx.txt 
usearch -cluster_aggd  mx.txt   -treeout tmp.tre -linkage min
tr -d "\n" < tmp.tre | sed "-es/zOTU_/'zOTU_/g" | sed "-es/:/':/g" | sed "-es/)':/):/g" > Results_${string4}/trees/zOTU.tree 
rm mx.txt 
rm tmp.tre

#>>>>>>move ZOTUs

mkdir Results_${string4}/zOTUs
cp zotus2.fa Results_${string4}/zOTUs/zOTUs.txt
mv zotus2.fa Results_${string4}/zOTUs/zOTUs.fa

mkdir Results_${string4}/taxonomy
mv EtsintaxR.txt Results_${string4}/taxonomy/sintax_EZtaxon.txt




