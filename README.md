# Fastq_2_zOTU-table workflow

This workflow contains the traditional 16S_rRNA-pipeline used in our analyses at UCPH-FOOD.

Before running the workflow, make sure to follow these steps:

1) Make sure you have administrator rights | UNIX environment
2) Download the free version of usearch (v11 is recommended) | https://www.drive5.com/usearch/download.html | make usearch executable by running: 
```
chmod 775 usearchv11***
```
4) 


```
reverse_file=(*R2*.fastq)
for file1 in *R1*.fastq
do

#>>>>>>>>adjusting label<<<<<<<<<<<<<<<<<<<<
tput setaf 4
string=${file1}
string2=$(sed "-es/_R[0-9].*\.fastq//; s/combined_//" <<< $string)
string3=$(sed "-es/-/_/g" <<< $string2)
echo "sampleID: ${string3}"

string4=$(sed "-es/S_[0-9].//g" <<< $string3)
echo ""
echo "project name: ${string4}"
echo ""
tput sgr0

#>>>>>>>>merging_ovelaping_pairs<<<<<<<<<<<<<<<<
usearch10 -fastq_mergepairs ${file1} -reverse ${reverse_file[i++]} -fastqout merged-${file1}.fq -fastq_minovlen 100 -relabel ${string4}.
#staggered reads arte removed by default in usearch10
done

tput setaf 4
echo ""
echo ">>>>merging fastq together<<<<"
#####################################cat merged-combined*fq > merged.fq
echo ""

#cleaning
#####################################rm merged-combined*.fq
mkdir ${string4}
mv *.fastq ${string4}/


#>>>>>>>Quality check<<<<<<<<<<<<<<<<<<<<<<<<<<<
echo ""
echo ">>>>checking fastq quality<<<<"
echo ""
tput sgr0
for file in merged*fq
do
usearch61 -fastq_filter $file -fastq_maxee 2.0  -fastq_truncqual 4 -fastq_minlen 130 -fastaout ${file}.fna
done
tput setaf 4
echo ""
echo ">>>>merged fastq filtered by quality<<<<"
tput sgr0

cat merged*.fna > qual-filt.fna
rm merged*.fna
rm merged*.fq

#>>>>>>>Find unique read sequences<<<<<<<<<<<<<<
tput setaf 4
echo ""
echo ">>>>searching unique sequences<<<<"
tput sgr0
usearch10 -fastx_uniques qual-filt.fna -sizeout -fastaout uniques.fa


usearch10 -unoise3 uniques.fa -zotus zotus.fa
sed "s/Zotu/zOTU_/g" zotus.fa > zotus2.fa 
rm zotus.fa

usearch10 -otutab qual-filt.fna -otus zotus2.fa -otutabout raw_zOTU.txt  -id 0.999
cat raw_zOTU.txt | sed 's/#OTU ID/OTUID/g' > raw_zOTU2.txt 
rm raw_zOTU.txt

usearch10 -sintax zotus2.fa -db '/mnt/ed058f9c-275f-4599-a283-e4e5927a7946/Databases_16s/rdp_16s_v16s_sp.fa' -strand both -tabbedout sintax0.txt -sintax_cutoff 0.8
usearch10 -sintax zotus2.fa -db '/mnt/ed058f9c-275f-4599-a283-e4e5927a7946/Databases_16s/Ez-taxon-format-sintax.fasta'  -strand both -tabbedout sintaxEt.txt -sintax_cutoff 0.8
TAB=$'\t'
cat sintax0.txt | awk '{print $1,$4}'| sed 's/"//g' | sed 's/p:/p__/g' | sed 's/c:/c__/g' | sed 's/o:/o__/g' | sed 's/f:/f__/g' | sed 's/g:/g__/g' | sed 's/s:/s__/g' | sed 's/,/;/g' | sed 's/ d:Bacteria;/'"${TAB}k__Bacteria;"'/g' | sed 's/ d:Archaea;/'"${TAB}k__Archaea;"'/g' > sintax.txt 

cat sintaxEt.txt | awk -F "\t" '{print $1,$4}'| sed 's/"//g' | sed 's/p:/p__/g' | sed 's/c:/c__/g' | sed 's/o:/o__/g' | sed 's/f:/f__/g' | sed 's/g:/g__/g' | sed 's/s:/s__/g' | sed 's/,/;/g' | sed 's/ d:Bacteria;/'"${TAB}k__Bacteria;"'/g' | sed 's/d:Bacteria/'"${TAB}k__Bacteria;"'/g' | sed 's/ d:Archaea;/'"${TAB}k__Archaea;"'/g' > Etsintax.txt 


rm sintax0.txt
rm sintaxEt.txt

echo "OTUID$TAB taxonomy" > header.txt 
cat header.txt | sed 's/ //g' > headers.txt 
cat headers.txt sintax.txt > sintaxR.txt
cat headers.txt Etsintax.txt > EtsintaxR.txt
rm headers.txt
rm header.txt
sed 's/#OTU ID/OTUID/g' < raw_zOTU2.txt > raw_zOTU_R2.txt
rm raw_zOTU2.txt


mkdir Results_${string4}
mkdir Results_${string4}/OTU-tables


echo '#!/usr/bin/env Rscript
X = read.table("raw_zOTU_R2.txt", header = TRUE)
Y <- read.table("sintaxR.txt", header = TRUE, fill = TRUE, na.string = "")
table <- merge(X,Y, by = "OTUID")
write.table(table, "OTU_completed.txt", sep="\t", row.names = F)' > biom.R

chmod 755 biom.R
./biom.R


echo '#!/usr/bin/env Rscript
X = read.table("raw_zOTU_R2.txt", header = TRUE)
Y <- read.table("EtsintaxR.txt", header = TRUE, fill = TRUE, na.string = "", sep="\t")
table <- merge(X,Y, by = "OTUID")
write.table(table, "OTU_completedEt.txt", sep="\t", row.names = F)' > biom3.R

chmod 755 biom3.R
./biom3.R


cat OTU_completed.txt | sed 's/"//g' | sed 's/OTUID/#OTU ID/g' | sed 's/NA/Unassigned/g' > zOTU_table_sintax.txt
biom convert -i zOTU_table_sintax.txt -o Results_${string4}/OTU-tables/zOTU_table_sintax.biom --process-obs-metadata taxonomy --table-type="OTU table" --to-json

cat OTU_completedEt.txt | sed 's/"//g' | sed 's/OTUID/#OTU ID/g' | sed 's/+/Unassigned/g' > zOTU_table_sintaxEt.txt
biom convert -i zOTU_table_sintaxEt.txt -o Results_${string4}/OTU-tables/zOTU_table_sintaxEt.biom --process-obs-metadata taxonomy --table-type="OTU table" --to-json

mv zOTU_table_sintax.txt Results_${string4}/OTU-tables/
mv zOTU_table_sintaxEt.txt Results_${string4}/OTU-tables/

tput setaf 4
echo ""
echo ">>>>generating phylogenetic tree<<<<"
tput sgr0
mkdir Results_${string4}/trees
usearch10 -cluster_agg zotus2.fa  -treeout tmp.tre -id 0.75 -linkage min


tr -d "\n" < tmp.tre | sed "-es/zOTU_/'zOTU_/g" | sed "-es/:/':/g" | sed "-es/)':/):/g" > Results_${string4}/trees/zOTU.tree 

mkdir Results_${string4}/zOTUs


assign_taxonomy.py -i zotus2.fa -o assign_output -r /mnt/ed058f9c-275f-4599-a283-e4e5927a7946/Databases_16s/gg_13_8_otus/rep_set/97_otus.fasta  -t /mnt/ed058f9c-275f-4599-a283-e4e5927a7946/Databases_16s/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt

echo "OTUID$TAB taxonomy$TAB hola$TAB tigre" > header.txt 
cat header.txt | sed 's/ //g' > headers.txt
cat headers.txt assign_output/zotus2_tax_assignments.txt | sed 's/ //g' > greengenesR.txt
rm headers.txt
rm header.txt

echo '#!/usr/bin/env Rscript
X = read.table("raw_zOTU_R2.txt", header = TRUE)
Z <- read.table("greengenesR.txt", header = TRUE, fill = TRUE, na.string = "")
Y = Z[,c(1,2)]
table <- merge(X,Y, by = "OTUID")
write.table(table, "OTU_completed.txt", sep="\t", row.names = F)' > biom.R

chmod 755 biom.R
./biom.R

cat OTU_completed.txt | sed 's/"//g' | sed 's/OTUID/#OTU ID/g' | sed 's/NA/Unassigned/g' > zOTU_table_GG.txt
biom convert -i zOTU_table_GG.txt -o Results_${string4}/OTU-tables/zOTU_table_GG.biom --process-obs-metadata taxonomy --table-type="OTU table" --to-json

#>>>>>>move ZOTUs
mv zOTU_table_GG.txt Results_${string4}/OTU-tables/
cp zotus2.fa Results_${string4}/zOTUs/zOTUs.txt
mv zotus2.fa Results_${string4}/zOTUs/zOTUs.fa

mkdir Results_${string4}/taxonomy
mv sintaxR.txt Results_${string4}/taxonomy/sintax_taxaRDP.txt
mv greengenesR.txt Results_${string4}/taxonomy/greengenes_taxa.txt
mv EtsintaxR.txt Results_${string4}/taxonomy/sintax_taxaEt.txt

rm -r assign_output/

notify-send "ANALYSIS COMPLEATED. Thank you! Krych & Castro-Mejía ;-) 2018"



tput setaf 4
echo ""
read -p ">>>>>>Would you like to remove half-products  (y/n) ?"  choice
if [ "$choice" == "y" ]; then


rm qual-filt.fna 
rm uniques.fa
rm sintax.txt
rm tmp.tre
rm biom*.R
rm OTU_completed.txt
rm raw_zOTU_R2.txt
rm Etsintax*
rm OTU_completedEt.txt

elif [ "$choice" == "n" ]; then
echo " O.K. files will be kept"
fi
tput sgr0
```
