#Script to obtain lncRNAs from mouse random transposn insertion screens after filtering for protein coding genes
#Adrienne Vancura (adrienne.vancura@dbmr.unibe.ch)

#variable 1: date, variable 2: setting liftover, variabe 3: setting liftover,  variable 4: which CIS to analyze (ALL)
#run from command line: sh CLC2:screenanalysis_NARCancer.sh 13092018 0.1 0.1 ALL 



projectname=$1
output_folder="$projectname.output"
mkdir /$output_folder$2$3$4$5
cd $output_folder$2$3$4$5


#download GENCODE genes from desired GENCODE version (https://www.gencodegenes.org) 
GENCODE="gencodeGenes.txt"

#download CCGD file from (http://ccgd-starrlab.oit.umn.edu/about.php)
#manually extract column T with CIS address, as command line tools don't extract the whole column
CCGD_file="CCGD_CISelements.txt"

#get mouse elements download from liftover instruction page (https://genome.ucsc.edu/cgi-bin/hgLiftOver)
MM10toHG19_file="mm10ToHg38.over.chain"
H19toMM10_file="hg38ToMm10.over.chain"



#get rid of empty lines, add chr to beginning of the line, make tabs instead of : and -
#subtract 1nt to get bed file for insertion of 1nt and no 0! with awk command and remove first row 

sed -e'/^\s*$/d; s/^/chr/; s///g; s/:/     /g; s/-/     /g' $CCGD_file | sed -e "1d" > CISelements.bed
echo "all CIS elements:"
cat CISelements.bed | wc -l
#22490


# #change input file so we can get input and output locations before and after liftover
# awk 'BEGIN{OFS="\t";} {print $1, $2, $3, $4 "input"  "\t" $1, $2, $3}' CISelements.bed > CISelementsALT.bed
# #make 0 based with $2-1 for proper liftover
# awk 'BEGIN{OFS="\t";} {print $1, ($2-1), $3, $4 ":" $5 ":" $6 ":" $7;}' CISelementsALT.bed > CISelementsALT1.bed


#change input file so we can get input and output locations before and after liftover
awk '{print $1, $2, $3, $4 "input" $1":"$2"-"$3}' CISelements.bed > CISelementsALT.bed
#awk '{print $1, $2, $3, $1, $2, $3}' CISelements.bed > CISelementsALT.bed

#make 0 based with $2-1 for proper liftover
awk '{print $1, ($2-1), $3, $4 }' CISelementsALT.bed > CISelementsALT1.bed
awk '{ print $1, $2, $3, $4, $3 - $2; }' CISelementsALT1.bed | awk '{print $1, $2, $3, $4}' | sed '/chrX, -1 inputchrX,:Y- 1/d' > CISelementsALT2.bed


ALLCIS_file="CISelementsALT2.bed"



####################################################################################
#liftOver
####################################################################################

#LiftOver default settings
#-minMatch=0.95 -minBlocks=0.1 for same specie
#-minMatch=0.1 -minBlocks=0.1 for inter specie

########################################################################################################################################################################

#liftOver -minMatch=$2 -minBlocks=$3
#to make it command line ready
#liftOver -minMatch=$2 -minBlocks=$3

variable4=$4

if [ "$4" == "ALL" ]
	then 
		liftOver -minMatch=$2 -minBlocks=$3 $ALLCIS_file $MM10toHG19_file outputALL.bed unliftedALL.bed
echo "ALL CIS with minMatch $2 and minBlocks $3"
cat outputALL.bed | wc -l
fi


####################################################################################
#bedtools merge to combind hCIS that are overlapping
####################################################################################
#input is 20466

#bedtools merge requires that you presort your data by chromosome and then by start position 

#$4 for ALL or EXACT or EXACTmod

sort -k1,1 -k2,2n output$4.bed > output.sorted.bed

#use bedtools merge to find overlapping hCIS and merge
bedtools merge -i output.sorted.bed > outputmerged.bed
echo "all sorted and merged hCIS $4:"
cat outputmerged.bed | wc -l
#6297



#prepare GENCODE bed file with all genes
#get the rows with gene in it, and not exon and transcript
#gencode from june 2018, check for newer versions too
GENCODEGTF="gencode.v28.annotation.gtf "

awk '$3 == "gene"' $GENCODEGTF > GENCODE.txt
echo "all GENCODE ids:"
cat $GENCODE | wc -l
#223830

#get rid of chrM
awk '!($1 == "chrM")' GENCODE.txt > GENCODE1.txt

awk '{print $1,$4,$5,$10}' GENCODE1.txt  > GENCODEids.bed

awk '{print $1,$4,$5}' GENCODE1.txt  > GENCODE.bed
echo "all GENCODE ids without chrM:"
cat GENCODEids.bed | wc -l
#58344


#remove 1 from start to make it 0 based as the other files
#they already should be O based as its a bed file
awk '{print $1, ($2-1), $3}' GENCODE.bed  > GENCODEm1.bed

#use tab with gsed, sed with tab character doesnt work otherwise
#only with genomic coordiantes
gsed -e's/ /\t/g' GENCODEm1.bed >GENCODEall.bed


#continue with this file as it has the genomic location and the ENSG
gsed -e's/ /\t/g' GENCODEids.bed >GENCODEallids.bed



####################################################################################
#intersect bed to get human ENSGs from genomic coordinates of hCIS elements
####################################################################################

#general input for intersect BED: intersectBed -a $MOUSE_FILE -b gencode.v24.intergenic.bed | sort -u | wc -l
#presort your data by chromosome and then by start position 

sort -k1,1 -k2,2n GENCODEallids.bed > GENCODEallids.sorted.bed

echo "input GENCODE and all hCIS:"
cat GENCODEallids.sorted.bed | wc -l
#58344
cat output$4.bed | wc -l
#20466



#keep input and output and add ENSG
intersectBed -wa -wb -a output$4.bed -b GENCODEallids.sorted.bed > overlapCISandGENCODE.bed 


echo "how many hCIS are found in GENCODE ids:"
awk '{print $7}' overlapCISandGENCODE.bed | sort | uniq -c | wc -l

echo "from a total of hCIS:"
cat output$4.bed | wc -l 


#compare to known functional cancer lncRNAs in the CLC2
CLC2GENES="CLC2018_ENSGonly.txt"


#get all ENSGs from liftover
awk '{print $8}' overlapCISandGENCODE.bed | cut -f1 -d"." | sed 's/"//g' | sort | uniq -c | awk '{print $2}' > ENSGsLiftOver.txt


echo "all ENSGs from liftover:"
cat ENSGsLiftOver.txt | wc -l
awk 'FNR==NR{a[$1];next}($1 in a){print}' $CLC2GENES  ENSGsLiftOver.txt > overlapCLC2andLiftoverResults.txt
echo "ENSGs which are in CLC2:"
cat overlapCLC2andLiftoverResults.txt | wc -l




####################################################################################
#address Gene types
####################################################################################


#prepare GENCODE bed file with all genes
#get the rows with gene in it, and not exon and transcript
GENCODEGTF="gencode.v28.annotation.gtf "

#setwd="/Users/ninavan/Dropbox/CLC2analysis/CISanalysis/20062018"


#only get lines with "gene" and get rid of chrM
awk '$3 == "gene"' $GENCODEGTF | awk '!($1 == "chrM")'> GENCODE.txt
echo "GENCODE" 
cat GENCODE.txt | wc -l
#58344

awk '{ print $1, $4, $5, $10 }' GENCODE.txt| sed 's/"//g' | sed 's/;//g' | cut -f1 -d"."  >  GENCODE.bed



####################################################################################
#get all protein_coding genes
####################################################################################

#we have to replace " and ; otherwise awk wont find protein:coding
cat GENCODE.txt | sed 's/"//g' | sed 's/;//g'> GENCODEnew.txt
awk '$12 == "protein_coding"' GENCODEnew.txt > proteinCodingGenes.txt
echo "protein coding genes"
cat proteinCodingGenes.txt | wc -l
#19888

#make bed file for all protein coding genes
awk '{ print $1, $4, $5, $10 }' proteinCodingGenes.txt >  proteinCodingGenesTEMP.bed

#subtract 1nt to get bed file for insertion of 1nt and no 0! with awk command and remove everything after the "."
gsed 's/ /\t/g' proteinCodingGenesTEMP.bed | awk '{print $1, ($2-1), $3, $4}' | cut -f1 -d"." > proteinCodingGenes.bed
#cat proteinCodingGenes.bed | wc -l
rm proteinCodingGenesTEMP.bed


####################################################################################
#get all CGC genes
####################################################################################

#download tsv from https://cancer.sanger.ac.uk/census

CENSUSALL="Census_all.tsv"
#get all ENSGs from this list
grep -o '\bENSG\w*' $CENSUSALL> CGCgenes.txt
#cat CGCgenes.txt | wc -l
#700

#get all CGC from proteincodingGenes.bed
grep -Ff CGCgenes.txt proteinCodingGenes.bed > CGCgenes.bed
echo "CGC genes in protein coding:" 
cat CGCgenes.bed | wc -l
#695

#some CGC genes are not protein coding
#go look for them in GENCODE file


grep -Ff CGCgenes.txt GENCODE.bed > CGCgenes1.bed
echo "CGC genes in GENCODE:" 
cat CGCgenes1.bed | wc -l
#696, only one more, whats with the others?

awk '{print $4}' CGCgenes1.bed > CGCgeneIDs.txt
echo "genes which are not in GENCODE, probably with new gene ID:"
awk 'FNR==NR{a[$1];next}!($1 in a){print}' CGCgeneIDs.txt CGCgenes.txt 
#ENSG00000258389 in v20 > v28 ENSG00000280757
#ENSG00000124693 is HIST1H3B > ENSG00000274267
#ENSG00000108292 is MLLT6 > ENSG00000275851
#ENSG00000140660 is RUNDC2A,SNX29  > ENSG00000048471 
#ENSG00000204645 is SSX4 > ENSG00000268009
#ENSG00000172660 is TAF15 > ENSG00000276833


####################################################################################
#get all nonCGC genes
####################################################################################


#to select only the ones not in both files
grep -vf CGCgenes.txt proteinCodingGenes.bed > nonCGCgenes.bed
echo "nonCGC genes"
cat nonCGCgenes.bed | wc -l
#19193


####################################################################################
#get all lnc genes
####################################################################################

#downlaod gtf from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.long_noncoding_RNAs.gtf.gz
#maybe get a newer version

GENCODElncRNA="gencode.v28.long_noncoding_RNAs.gtf"

#get only row with ENSGs, get rid of " and ;, get rid of everything after the "." and only get uniq ones 
awk '{print $10}' $GENCODElncRNA | sed 's/"//g' | sed 's/;//g' | cut -f1 -d"." | sort | uniq -c | awk '{print $2}' | sed -e "1d"| sed -e "1d" > lncRNAgenes.txt
#cat lncRNAgenes.txt | wc -l
#15767


#get all lncRNA from GENCODE.txt
grep -Ff lncRNAgenes.txt GENCODE.bed > lncRNA.bed
#cat lncRNA.bed | awk '{print $4}'| sort | uniq -c | wc -l
gsed -e 's/ /\t/g' lncRNA.bed > lncRNA1.bed


####################################################################################
#get all CLC genes
####################################################################################

#new CLC gene list (literature lncRNAs from this version)

CLC2genes="CLC2018_ENSGonly.txt"
#cat CLC2genes.txt | wc -l
#244

#get all CLC from lncRNA.bed
grep -Ff $CLC2genes lncRNA.bed > CLC2genes.bed
echo "CLC2 genes in lncRNA annotation file:" 
cat CLC2genes.bed | wc -l
#208
#some CLC2 genes are not in the lncRNA annotation file, check in the gencode file


grep -Ff $CLC2genes GENCODE.bed > CLC2genes1.bed
echo "CLC2 genes in GENCODE:" 
cat CLC2genes1.bed | wc -l
#243
awk '{print $4}' CLC2genes1.bed > CLC2genes1.txt
#cat CLC2genes1.txt | wc -l



####################################################################################
#get all nonCLC genes
####################################################################################


#to select only the ones not in both files
#use lncRNA annotation file 
grep -vf $CLC2genes lncRNA.bed > nonCLC2genes.bed
echo "nonCLC genes"
cat nonCLC2genes.bed | wc -l
#15571



####################################################################################
#get whole genome
####################################################################################
#go to UCSC table browser and put in ALL tables for group and then ChromInfo as table to get genomic coordinates
#https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=677420819_Ap82DTAgagxrVxh9tvGk5OU3EiFc&clade=mammal&org=Human&db=hg38&hgta_group=allTables&hgta_track=hg38&hgta_table=chromInfo&hgta_regionType=genome&position=chr17_KI270857v1_alt%3A43987-81774&hgta_outputType=primaryTable&hgta_outFileName=


#cat wholeGenome.bed | wc -l
#23



####################################################################################
#delete overlapping regions 
####################################################################################

#important to use the bed files originating from protein coding or lncRNA bed files, otherwise we cant be sure about the origin

wholeGenome="/Users/ninavan/Dropbox/CLC2analysis/databases/wholeGenome/wholeGenome1.bed"
#wholeGenome=wholeGenome.bed
#CGC=CGCgenes.bed
#nonCGC=nonCGCgenes.bed
#proteinCoding=proteinCodingGenes.bed
#CLC2=CLC2genes1.bed
#nonCLC2=nonCLC2genes.bed
#lncRNA=lncRNA.bed

#make tab delimited file and sort
gsed -e's/ /\t/g' nonCGCgenes.bed | awk '{print $1,$2,$3}' | gsed -e's/ /\t/g' > nonCGCgenes1.bed
gsed -e's/ /\t/g' CGCgenes.bed | awk '{print $1,$2,$3}' | gsed -e's/ /\t/g' > CGCgenes1.bed
gsed -e's/ /\t/g' CLC2genes1.bed | awk '{print $1,$2,$3}' | gsed -e's/ /\t/g' > CLCgenes2.bed
gsed -e's/ /\t/g' nonCLC2genes.bed | awk '{print $1,$2,$3}' | gsed -e's/ /\t/g'> nonCLCgenes1.bed
gsed -e's/ /\t/g' $wholeGenome | awk '{print $1,$2,$3}' | gsed -e's/ /\t/g'> wholeGenome1.bed


#merge all overlapping elements in a file with bedtools merge
bedtools merge -i nonCGCgenes1.bed > nonCGCgenes1_M.bed
bedtools merge -i CGCgenes1.bed > CGCgenes1_M.bed
bedtools merge -i CLCgenes2.bed > CLCgenes2_M.bed
bedtools merge -i nonCLCgenes1.bed > nonCLCgenes1_M.bed
bedtools merge -i wholeGenome1.bed > wholeGenome1_M.bed



#hierarchy like this: CGC > nonCGC > CLC > nonCLC > intergenic.
#intersect the CIS elements to protein coding before noncoding section!


#1.nonCGC - CGC = correct-nonCGC
subtractBed -a nonCGCgenes1_M.bed -b CGCgenes1_M.bed > correct_nonCGC.bed


#2.CLC - CGC - nonCGC = correct-CLC
subtractBed -a CLCgenes2_M.bed -b CGCgenes1_M.bed > correct_CLC.bed
subtractBed -a correct_CLC.bed -b correct_nonCGC.bed > correct_CLC1.bed
rm correct_CLC.bed


#3.nonCLC - CLC - CGC - nonCGC = correct-nonCLC
subtractBed -a nonCLCgenes1_M.bed -b correct_CLC1.bed > correct_nonCLC.bed
subtractBed -a correct_nonCLC.bed -b CGCgenes1_M.bed > correct_nonCLC1.bed
subtractBed -a correct_nonCLC1.bed -b correct_nonCGC.bed > correct_nonCLC2.bed

#4.whole-genome - nonCLC - CLC - CGC - nonCGC = intergenic
subtractBed -a wholeGenome1_M.bed -b correct_nonCLC2.bed > correct_wholeGenome.bed
subtractBed -a correct_wholeGenome.bed -b correct_CLC1.bed > correct_wholeGenome1.bed
subtractBed -a correct_wholeGenome1.bed -b correct_nonCGC.bed > correct_wholeGenome2.bed
subtractBed -a correct_wholeGenome2.bed -b CGCgenes1_M.bed > correct_wholeGenome3.bed


####################################################################################
#end product
####################################################################################


#Now you will have non-ovelapping BED files, one for each of the five categories: CGC, nonCLC, CLC, nonCLC and intergenic.
#correct_nonCLC2.bed ,correct_CLC1.bed ,CGCgenes1.bed ,correct_nonCGC.bed,correct_wholeGenome3.bed

#Make sure that the sum of the length of each element in those five files adds up to the total size of the genome (which you can calculate from the “whole-genome” BED file).
awk '{SUM += $3-$2} END {print SUM}' correct_nonCLC2.bed > SUMbedFiles.txt
awk '{SUM += $3-$2} END {print SUM}' correct_CLC1.bed >> SUMbedFiles.txt
awk '{SUM += $3-$2} END {print SUM}' CGCgenes1_M.bed >> SUMbedFiles.txt
awk '{SUM += $3-$2} END {print SUM}' correct_nonCGC.bed >> SUMbedFiles.txt
awk '{SUM += $3-$2} END {print SUM}' correct_wholeGenome3.bed >> SUMbedFiles.txt
awk '{SUM += $1} END {print SUM}' SUMbedFiles.txt >> SUMbedFiles.txt
echo "SUM of all locations:"
cat  SUMbedFiles.txt
#3088269832
echo "SUM of whole genome:"
awk '{SUM += $3-$2} END {print SUM}' wholeGenome1.bed
#3088269832



####################################################################################
#prepare data for figure depicting how many CIS elements are in each category
####################################################################################

#use the hCIS BED file


#CISelements="/Users/ninavan/Dropbox/CLC2analysis/CISanalysis/14062018/outputALL.bed"

CISelements=output$4.bed


echo "nonCLC" > CISelementsGeneFamilies.txt
intersectBed -a $CISelements -b correct_nonCLC2.bed | sort -u | wc -l >> CISelementsGeneFamilies.txt
echo "CLC2" >> CISelementsGeneFamilies.txt
intersectBed -a $CISelements -b correct_CLC1.bed | sort -u | wc -l >> CISelementsGeneFamilies.txt
echo "CGC" >> CISelementsGeneFamilies.txt
intersectBed -a $CISelements -b CGCgenes1_M.bed | sort -u | wc -l >> CISelementsGeneFamilies.txt
echo "nonCGC" >> CISelementsGeneFamilies.txt
intersectBed -a $CISelements -b correct_nonCGC.bed | sort -u | wc -l >> CISelementsGeneFamilies.txt
echo "intergenic" >> CISelementsGeneFamilies.txt
intersectBed -a $CISelements -b correct_wholeGenome3.bed | sort -u | wc -l >> CISelementsGeneFamilies.txt


#some CIS can be mapped to two or more genes
#make file to a table
awk 'ORS=NR%2?" ":"\n"' CISelementsGeneFamilies.txt > CISelementsGeneFamilies1.txt

#sum of all hCIS in file (for percentage calculations)
#for details check: https://unix.stackexchange.com/questions/174371/calculate-and-divide-by-total-with-awk
awk 'FNR==NR{s+=$2;next;} {printf "%s\t%s\t%s%%\n",$1,$2,100*$2/s}' CISelementsGeneFamilies1.txt CISelementsGeneFamilies1.txt > CISelementsGeneFamiliesPercentages.txt



awk '{{SUM += $2}/100 * $2 END {print}}' CISelementsGeneFamilies1.txt



echo "CISelements in gene families:"
cat CISelementsGeneFamilies1.txt


#Finally, some hCIS may overlap two genes at the same type, like a CGC and a intergenic region. So you need to find out these cases. With a commands like these:

####################################################################################
#get exact hCIS
####################################################################################

#By default, bedtools intersect will report an overlap between A and B so long as there is at least one base pair is overlapping.
#so if 1 bp overlaps its added to the bin

#add extra code to sort out CIS elements that mapped to CGC
intersectBed -a $CISelements -b CGCgenes1_M.bed > CISelements_CGC.bed
intersectBed -v -a $CISelements -b CISelements_CGC.bed > CISelements_afterCGC.bed 
#intersectBed -v -a outputCISVAR.bed -b CISelements_CGC.bed > CISelements_afterCGC.bed 


intersectBed -a CISelements_afterCGC.bed -b correct_nonCGC.bed > CISelements_nonCGC.bed
intersectBed -v -a CISelements_afterCGC.bed -b CISelements_nonCGC.bed > CISelements_afterNonCGC.bed 


intersectBed -a CISelements_afterNonCGC.bed -b correct_CLC1.bed > CISelements_CLC2.bed
intersectBed -v -a CISelements_afterNonCGC.bed -b CISelements_CLC2.bed > CISelements_afterCLC2.bed 


intersectBed -a CISelements_afterCLC2.bed -b correct_nonCLC2.bed > CISelements_nonCLC2.bed
intersectBed -v -a CISelements_afterCLC2.bed -b CISelements_nonCLC2.bed > CISelements_afterNonCLC2.bed 


intersectBed -a CISelements_afterNonCLC2.bed -b correct_wholeGenome3.bed > CISelements_intergenic.bed
intersectBed -v -a CISelements_afterNonCLC2.bed -b CISelements_intergenic.bed > CISelements_afterIntergenic.bed 
cat CISelements_afterIntergenic.bed | wc -l
#at the end it should be 0


#input file: 20466, cat outputALL.bed | wc -l
#total 22898
#CGC hCIS: 2611
#nonCGC: 19801
#CLC2: 60
#nonCLC2: 197
#intergenic: 229



awk '{print $4}' CISelements_CGC.bed | sort | uniq -c | awk '{print $2}' | wc -l
awk '{print $4}' CISelements_nonCGC.bed | sort | uniq -c | awk '{print $2}' | wc -l
awk '{print $4}' CISelements_CLC2.bed | sort | uniq -c | awk '{print $2}' | wc -l
awk '{print $4}' CISelements_nonCLC2.bed | sort | uniq -c | awk '{print $2}' | wc -l
awk '{print $4}' CISelements_intergenic.bed | sort | uniq -c | awk '{print $2}' | wc -l


awk '{print $4}' outputALL.bed  | sort | uniq -c | awk '{print $2}' | wc -l
#16430




#unique
#input file: 20466, cat outputALL.bed | wc -l
#total 22898
#CGC hCIS: 1915
#nonCGC: 14112
#CLC2: 50
#nonCLC2: 144
#intergenic: 209
#unique hCIS: 16'430




awk '{print $4}' CISelements_CGC.bed | sort | uniq -c | awk '{print $2}' | gsed -e 's/:/\t/g' | gsed -e 's/-/\t/g' | awk '{print $3-$2}' > CISlengths_unique.txt
awk '{print $4}' CISelements_nonCGC.bed | sort | uniq -c | awk '{print $2}' | gsed -e 's/:/\t/g' | gsed -e 's/-/\t/g' | awk '{print $3-$2}' >> CISlengths_unique.txt


awk '{print $4}' CISelements_CLC2.bed | sort | uniq -c | awk '{print $2}' | gsed -e 's/:/\t/g' | gsed -e 's/-/\t/g' | awk '{print $3-$2}' >> CISlengths_unique.txt
awk '{print $4}' CISelements_nonCLC2.bed | sort | uniq -c | awk '{print $2}' | gsed -e 's/:/\t/g' | gsed -e 's/-/\t/g' | awk '{print $3-$2}' >> CISlengths_unique.txt
awk '{print $4}' CISelements_intergenic.bed | sort | uniq -c | awk '{print $2}' | gsed -e 's/:/\t/g' | gsed -e 's/-/\t/g' | awk '{print $3-$2}' >> CISlengths_unique.txt

cat CISlengths_unique.txt | sort | uniq -c 

#exact:
#1196 0

#16430 - 1196 = 15234


intersectBed -wa -wb -a CISelements_nonCLC2.bed -b lncRNA1.bed > CISelements_nonCLC2_ENSGs_lncRNA.bed 
intersectBed -wa -wb -a CISelements_CLC2.bed -b lncRNA1.bed > CISelements_CLC2_ENSGs_lncRNA.bed 



intersectBed -wa -wb -a CISelements_nonCLC2.bed -b GENCODEallids.bed > CISelements_nonCLC2_ENSGs.bed 
intersectBed -wa -wb -a CISelements_CLC2.bed -b GENCODEallids.bed > CISelements_CLC2_ENSGs.bed 
intersectBed -wa -wb -a CISelements_CGC.bed -b GENCODEallids.bed > CISelements_CGC_ENSGs.bed 
intersectBed -wa -wb -a CISelements_nonCGC.bed -b GENCODEallids.bed > CISelements_nonCGC_ENSGs.bed 
intersectBed -wa -wb -a CISelements_intergenic.bed -b GENCODEallids.bed > CISelements_intergenic_ENSGs.bed 



awk '{print $8}' CISelements_nonCLC2_ENSGs.bed | cut -f1 -d"." | sed 's/"//g' | sort | uniq -c | awk '{print $2}' > CISelements_nonCLC2_ENSGs.txt
awk '{print $8}' CISelements_CLC2_ENSGs.bed | cut -f1 -d"." | sed 's/"//g' | sort | uniq -c | awk '{print $2}' > CISelements_CLC2_ENSGs.txt
awk '{print $8}' CISelements_CGC_ENSGs.bed | cut -f1 -d"." | sed 's/"//g' | sort | uniq -c | awk '{print $2}' > CISelements_CGC_ENSGs.txt
awk '{print $8}' CISelements_nonCGC_ENSGs.bed | cut -f1 -d"." | sed 's/"//g' | sort | uniq -c | awk '{print $2}' > CISelements_nonCGC_ENSGs.txt
awk '{print $8}' CISelements_intergenic_ENSGs.bed | cut -f1 -d"." | sed 's/"//g' | sort | uniq -c | awk '{print $2}' > CISelements_intergenic_ENSGs.txt


awk '{print $8}' CISelements_CLC2_ENSGs_lncRNA.bed | cut -f1 -d"." | sed 's/"//g' | sort | uniq -c | awk '{print $2}' > CISelements_CLC2_ENSGs_lncRNAannotation.txt
awk '{print $8}' CISelements_nonCLC2_ENSGs_lncRNA.bed | cut -f1 -d"." | sed 's/"//g' | sort | uniq -c | awk '{print $2}' > CISelements_nonCLC2_ENSGs_lncRNAannotation.txt

 



echo "overlap of hCIS elements:"
comm -12 CISelements_nonCLC2.bed CISelements_CLC2.bed | wc -l 
comm -12 CISelements_nonCLC2.bed CISelements_CGC.bed | wc -l
comm -12 CISelements_nonCLC2.bed CISelements_nonCGC.bed | wc -l
comm -12 CISelements_nonCLC2.bed CISelements_intergenic.bed | wc -l

comm -12 CISelements_CLC2.bed CISelements_CGC.bed | wc -l 
comm -12 CISelements_CLC2.bed CISelements_nonCGC.bed | wc -l 
comm -12 CISelements_CLC2.bed CISelements_intergenic.bed | wc -l 

comm -12 CISelements_CGC.bed CISelements_nonCGC.bed | wc -l 
comm -12 CISelements_CGC.bed CISelements_intergenic.bed | wc -l 

comm -12 CISelements_nonCGC.bed CISelements_intergenic.bed | wc -l 

#all results in 0


########################################


#make file with exact CIS from mouse $4, exact location in human $5,$6,$7, also calculate size of hCIS with $7-$6, and ENSG $8


CCGDselection="/Users/ninavan/Dropbox/CLC2analysis/CISanalysis/CCGDselection.txt"

#change format of CCGD file so it works with AWK
tr "\r" "\n" < $CCGDselection > CCGDselection1.txt

sed 's/ //g' CCGDselection1.txt | awk 'NF==4' | sort > CCGDselection2.txt

cat CCGDselection2.txt | wc -l
#all mouse CIS elements, 22491, correct


#remove input:chr so columns match with CCGD selection

#cat overlapCISandGENCODE.bed | sed 's/inputchr//g' | gsed -e's/ /:/g' | sed 's/:/-/2' | sed 's/"//g' | sed 's/;//g' | cut -f1 -d"." > overlapCISandGENCODE.txt

#awk '{print $4, $1, $2, $3, $5, $6, $7, $8;}' overlapCISandGENCODE.txt | sort > overlapCISandGENCODE1.txt


#change this section, as the overlapCISandGENCODE1.txt file contains several ENSGs per CIS, messing up protein coding ones and non coding 
#join CCGDselection2.txt overlapCISandGENCODE1.txt > CISelementsALLINFO.txt


cat CISelements_nonCLC2_ENSGs.bed | sed 's/inputchr//g' | gsed -e's/ /:/g' | sed 's/:/-/2' | sed 's/"//g' | sed 's/;//g' | cut -f1 -d"." > CISelements_nonCLC2_ENSGsModified.txt
awk '{print $4, $1, $2, $3, $5, $6, $7, $8;}' CISelements_nonCLC2_ENSGsModified.txt | sort > CISelements_nonCLC2_ENSGsModified1.txt

cat CISelements_CLC2_ENSGs.bed | sed 's/inputchr//g' | gsed -e's/ /:/g' | sed 's/:/-/2' | sed 's/"//g' | sed 's/;//g' | cut -f1 -d"." > CISelements_CLC2_ENSGsModified.txt
awk '{print $4, $1, $2, $3, $5, $6, $7, $8;}' CISelements_CLC2_ENSGsModified.txt | sort > CISelements_CLC2_ENSGsModified1.txt




join CCGDselection2.txt CISelements_nonCLC2_ENSGsModified1.txt > CISelements_nonCLC2_ALLINFO.txt
join CCGDselection2.txt CISelements_CLC2_ENSGsModified1.txt > CISelements_CLC2_ALLINFO.txt


# mouse CIS, cancer type, gain loss, human CIS chr start end size 	gene ID gene name   paste variable $4 (EXACT,mod, varlen)
#awk 'BEGIN{OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, ($7-$6) }' CISelementsALLINFO.txt > CISelementsALLINFO_1.txt
# awk  -v var="exact" '$12 <= var' CISelementsALLINFO.txt | awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' CISelementsALLINFO.txt

awk 'BEGIN{OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, ($7-$6) }' CISelements_nonCLC2_ALLINFO.txt > CISelements_nonCLC2_ALLINFO_1.txt
awk 'BEGIN{OFS="\t";} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, ($7-$6) }' CISelements_CLC2_ALLINFO.txt > CISelements_CLC2_ALLINFO_1.txt




#only get CLC and nonCLC from that file
#be careful! the CISelementsALLINFO is a "dirty" file with overlaying protein and non coding genes! 


#grep -Ff CISelements_CLC2_ENSGs.txt CISelementsALLINFO_1.txt | sort > CISelementsALLINFO_CLC.txt
#grep -Ff CISelements_nonCLC2_ENSGs.txt CISelementsALLINFO_1.txt |sort > CISelementsALLINFO_nonCLC.txt

grep -Ff CISelements_CLC2_ENSGs_lncRNAannotation.txt CISelements_CLC2_ALLINFO_1.txt | sort > CISelementsALLINFO_CLC_lncRNAannotation.txt
grep -Ff CISelements_nonCLC2_ENSGs_lncRNAannotation.txt CISelements_nonCLC2_ALLINFO_1.txt | sort > CISelementsALLINFO_nonCLC_lncRNAannotation.txt




#awk '{print $11, $5, $6, $7, $8, $9, $10, $12, $1, $2, $3, $4;}' CISelementsALLINFO_CLC.txt | sort > CISelementsALLINFO_CLC_new.txt
#awk '{print $11, $5, $6, $7, $8, $9, $10, $12, $1, $2, $3, $4;}' CISelementsALLINFO_nonCLC.txt | sort > CISelementsALLINFO_nonCLC_new.txt

awk '{print $11, $5, $6, $7, $8, $9, $10, $12, $1, $2, $3, $4;}' CISelementsALLINFO_nonCLC_lncRNAannotation.txt | sort > CISelementsALLINFO_nonCLC_new_lncRNAannotation.txt

awk '{print $11, $5, $6, $7, $8, $9, $10, $12, $1, $2, $3, $4;}' CISelementsALLINFO_CLC_lncRNAannotation.txt | sort > CISelementsALLINFO_CLC_new_lncRNAannotation.txt

#run Rscript to get Ensemble gene names 
Rscript /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/scripts/biomart.R


sed 's/"//g' geneNames_CLC.txt | awk '{print $2, $3;}' | sort > geneNames_CLC1.txt
sed 's/"//g' geneNames_nonCLC.txt | awk '{print $2, $3;}'| sort > geneNames_nonCLC1.txt

 
#join geneNames_CLC1.txt CISelementsALLINFO_CLC_new.txt > CISelementsALLINFO_CLC_FINAL.txt
#join geneNames_nonCLC1.txt CISelementsALLINFO_nonCLC_new.txt > CISelementsALLINFO_nonCLC_FINAL.txt


join geneNames_nonCLC1.txt CISelementsALLINFO_nonCLC_new_lncRNAannotation.txt > CISelementsALLINFO_nonCLC_FINAL_lncRNAannotation.txt
join geneNames_CLC1.txt CISelementsALLINFO_CLC_new_lncRNAannotation.txt > CISelementsALLINFO_CLC_FINAL_lncRNAannotation.txt



#make bed file with CIS locations to check in genome browser
#delete doubles one with the sort function

awk '{print chr$3,$4,$5}' CISelementsALLINFO_nonCLC_FINAL_lncRNAannotation.txt | sort -u -k1,1 -k2,2 -k3,3 -s > hCISforGenomeBrowser_nonCLC.bed

awk '{print chr$3,$4,$5}' CISelementsALLINFO_CLC_FINAL_lncRNAannotation.txt | sort -u -k1,1 -k2,2 -k3,3 -s > hCISforGenomeBrowser_CLC.bed




#make one file for genome browser with all CIS elements to be able to check quickly!
awk '{print chr$1,$2,$3}' CISelements_nonCLC2.bed | sort -u -k1,1 -k2,2 -k3,3 -s > hCISforGenomeBrowser_nonCLC_all.bed
awk '{print chr$1,$2,$3}' CISelements_CLC2.bed | sort -u -k1,1 -k2,2 -k3,3 -s > hCISforGenomeBrowser_CLC_all.bed
awk '{print chr$1,$2,$3}' CISelements_CGC.bed | sort -u -k1,1 -k2,2 -k3,3 -s > hCISforGenomeBrowser_CGC_all.bed
awk '{print chr$1,$2,$3}' CISelements_nonCGC.bed | sort -u -k1,1 -k2,2 -k3,3 -s > hCISforGenomeBrowser_nonCGC_all.bed
awk '{print chr$1,$2,$3}' CISelements_intergenic.bed | sort -u -k1,1 -k2,2 -k3,3 -s > hCISforGenomeBrowser_intergenic_all.bed


cat hCISforGenomeBrowser_nonCLC_all.bed hCISforGenomeBrowser_CLC_all.bed hCISforGenomeBrowser_CGC_all.bed hCISforGenomeBrowser_nonCGC_all.bed hCISforGenomeBrowser_intergenic_all.bed > hCISforGenomeBrowser_all.bed


cat hCISforGenomeBrowser_all.bed | gsed -e's/ /\t/g' > hCISforGenomeBrowser_all_formatted.bed


#cat hCISforGenomeBrowser_nonCLC.bed | sort -u -k1,1 -k2,2 -k3,3 -s 





#######################################################################

###################### M O U S E - S E C T I O N ######################

#######################################################################

#check which CIS elements in mouse map to protein coding mouse genes
#check which genes are at CIS elements from mouse







####################################################################################

GENCODEMOUSE="/Users/ninavan/Dropbox/CLC2analysis/databases/gencode/gencode.vM19.annotation.gtf"


#do the same for the mouse GENCODE
awk '$3 == "gene"' $GENCODEMOUSE | awk '!($1 == "chrM")'> GENCODEmouse.txt
echo "GENCODE mouse:" 
cat GENCODEmouse.txt | wc -l
#54409

awk '{ print $1, $4, $5, $10 }' GENCODEmouse.txt| sed 's/"//g' | sed 's/;//g' | cut -f1 -d"."  >  GENCODEmouse.bed


####################################################################################
#get all protein_coding MOUSE genes
####################################################################################

#we have to replace " and ; otherwise awk wont find protein:coding
cat GENCODEmouse.txt | sed 's/"//g' | sed 's/;//g'> GENCODEmousenew.txt
awk '$12 == "protein_coding"' GENCODEmousenew.txt > proteinCodingGenesMouse.txt
echo "protein coding genes"
cat proteinCodingGenesMouse.txt | wc -l
#21956

#make bed file for all protein coding genes
awk '{ print $1, $4, $5, $10 }' proteinCodingGenesMouse.txt > proteinCodingGenesTEMPMouse.bed

#subtract 1nt to get bed file for insertion of 1nt and no 0! with awk command and remove everything after the "."
gsed 's/ /\t/g' proteinCodingGenesTEMPMouse.bed | awk '{print $1, ($2-1), $3, $4}' | cut -f1 -d"." | gsed 's/ /\t/g' > proteinCodingGenesMouse.bed
#cat proteinCodingGenes.bed | wc -l
rm proteinCodingGenesTEMPMouse.bed


awk '$3 == "gene"' GENCODEmousenew.txt > allMouseGenes.txt
gsed 's/ /\t/g' allMouseGenes.txt | awk '{print $1, ($4-1), $5, $10}' | cut -f1 -d"." | gsed 's/ /\t/g' > GencodeGenesMouseALL.bed







####################################################################################
#get all lnc genes MOUSE to be able to sort out CIS in CLC and nonCLC but mapping to protein coding in mouse!!!
####################################################################################




#downlaod gtf from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.long_noncoding_RNAs.gtf.gz
GENCODElncRNA="/Users/ninavan/Dropbox/CLC2analysis/CISanalysis/20062018/gencode.v28.long_noncoding_RNAs.gtf"
GENCODElncRNAMOUSE="/Users/ninavan/Dropbox/CLC2analysis/databases/gencode/gencode.vM19.long_noncoding_RNAs.gtf"


#get only row with ENSGs, get rid of " and ;, get rid of everything after the "." and only get uniq ones 
awk '{print $10}' $GENCODElncRNAMOUSE | sed 's/"//g' | sed 's/;//g' | cut -f1 -d"." | sort | uniq -c | awk '{print $2}' | sed -e "1d"| sed -e "1d"  > lncRNAgenesMOUSE.txt
#cat lncRNAgenes.txt | wc -l
#15767

#get all lncRNA from GENCODE.txt
grep -Ff lncRNAgenesMOUSE.txt GENCODEmouse.bed | gsed 's/ /\t/g' > lncRNAMouse.bed
#cat lncRNA.bed | awk '{print $4}'| sort | uniq -c | wc -l
#some of the IDs have two genomic locations? why?
echo "lncRNA mouse:" 

cat lncRNAMouse.bed | wc -l
#12839




#get CLC and nonCLC CIS hits

awk '{print $9}'  CISelementsALLINFO_nonCLC_new_lncRNAannotation.txt | sort | uniq -c | awk '{print $2}' | sed -e' s/^/chr/; s/:/     /g; s/-/     /g' | sort > nonCLC_mouseCIS.bed
#these are 150 in total 
awk '{print $9}'  CISelementsALLINFO_CLC_new_lncRNAannotation.txt | sort | uniq -c | awk '{print $2}' | sed -e' s/^/chr/; s/:/     /g; s/-/     /g' | sort > CLC_mouseCIS.bed



cat nonCLC_mouseCIS.bed | gsed 's/     /\t/g' > nonCLC_mouseCIS_corr.bed
cat CLC_mouseCIS.bed | gsed 's/     /\t/g' > CLC_mouseCIS_corr.bed

cat proteinCodingGenesMouse.bed | gsed 's/ /\t/g' > proteinCodingGenesMouse_corr.bed



intersectBed -wa -wb -a nonCLC_mouseCIS_corr.bed -b proteinCodingGenesMouse_corr.bed > CISelements_nonCLC_inProteinCoding_MOUSE.bed

intersectBed -wa -wb -a CLC_mouseCIS_corr.bed -b proteinCodingGenesMouse_corr.bed > CISelements_CLC_inProteinCoding_MOUSE.bed


awk '{print $1":"$2"-"$3}' CISelements_nonCLC_inProteinCoding_MOUSE.bed | sed 's/chr//g;' > CISelements_nonCLC_inProteinCoding_MOUSE_CISToRemove.txt
awk '{print $1":"$2"-"$3}' CISelements_CLC_inProteinCoding_MOUSE.bed | sed 's/chr//g;' > CISelements_CLC_inProteinCoding_MOUSE_CISToRemove.txt
#44 and 3 CIS to remove 

#sort out these mouse CIS from the files
#first but CIS at the beginning for the code to work

awk '{print $9,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' CISelementsALLINFO_nonCLC_new_lncRNAannotation.txt | sort > CISelementsALLINFO_nonCLC_new_lncRNAannotation2.txt  
awk '{print $9,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' CISelementsALLINFO_CLC_new_lncRNAannotation.txt | sort > CISelementsALLINFO_CLC_new_lncRNAannotation2.txt  






grep -vf CISelements_nonCLC_inProteinCoding_MOUSE_CISToRemove.txt CISelementsALLINFO_nonCLC_new_lncRNAannotation2.txt | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > CISelementsALLINFO_nonCLC_new_lncRNAannotationFINAL.txt

grep -vf CISelements_CLC_inProteinCoding_MOUSE_CISToRemove.txt CISelementsALLINFO_CLC_new_lncRNAannotation2.txt | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > CISelementsALLINFO_CLC_new_lncRNAannotationFINAL.txt
















#only get unique ones, no duplicates
cat CISelementsALLINFO_nonCLC_new_lncRNAannotationFINAL.txt | sort | uniq -c | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > uniqueCISelements_nonCLC_lncRNAannotation_afterMousePC_FINAL.txt

cat uniqueCISelements_nonCLC_lncRNAannotation_afterMousePC_FINAL.txt | wc -l
#183

awk '{print $1}' uniqueCISelements_nonCLC_lncRNAannotation_afterMousePC_FINAL.txt | sort | uniq -c | wc -l
#get exact number of ENSGs
#107



cat CISelementsALLINFO_CLC_new_lncRNAannotationFINAL.txt | sort | uniq -c | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > uniqueCISelements_CLC_lncRNAannotation_afterMousePC_FINAL.txt

cat uniqueCISelements_CLC_lncRNAannotation_afterMousePC_FINAL.txt | wc -l
#78
awk '{print $1}' uniqueCISelements_CLC_lncRNAannotation_afterMousePC_FINAL.txt | sort | uniq -c | wc -l
#35











#do the same for the intergenic CIS 

awk '{print $4}' CISelements_intergenicSorted.bed | sort | uniq -c | awk '{print $2}' | sed -e'  s/input//g; s/:/     /g; s/-/     /g' | sort > intergenicHITS_mouseCIS.bed



#CISelements_intergenicSorted.bed 








##########################################


#Run the Rscript to generate barplots
Rscript /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/scripts/genefamilies.R


echo "getting exact hCIS"



####################################################################################
#copy output data in to new folder to check results more easily
####################################################################################

mkdir /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results

cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/GeneTypes.pdf /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/GeneTypes.pdf

cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/CISelementsGeneFamiliesPercentages.txt /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/CISelementsGeneFamiliesPercentages.txt

cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/CISelements_nonCLC2_ENSGs.txt /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/CISelements_nonCLC2_ENSGs.txt
cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/CISelements_CLC2_ENSGs.txt /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/CISelements_CLC2_ENSGs.txt 
cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/CISelements_CGC_ENSGs.txt  /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/CISelements_CGC_ENSGs.txt 
cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/CISelements_nonCGC_ENSGs.txt /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/CISelements_nonCGC_ENSGs.txt 
cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/CISelements_intergenic_ENSGs.txt  /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/CISelements_intergenic_ENSGs.txt 

cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/CISelementsALLINFO_CLC_FINAL.txt /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/CISelementsALLINFO_CLC_FINAL.txt
cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/CISelementsALLINFO_nonCLC_FINAL.txt  /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/CISelementsALLINFO_nonCLC_FINAL.txt




####################################################################################
#check if ENSGs overlap with PCAWG or EXINATOR genes
####################################################################################


PCAWGgenes="/Users/ninavan/Dropbox/CLC2analysis/databases/ExinatorPCAWG/PCAWG.txt"
Exinator1genes="/Users/ninavan/Dropbox/CLC2analysis/databases/ExinatorPCAWG/ExInAtor1.txt"
Exinator2genes="/Users/ninavan/Dropbox/CLC2analysis/databases/ExinatorPCAWG/ExInAtor2.txt"

grep -Ff $PCAWGgenes CISelements_nonCLC2_ENSGs.txt > overlapPCAWGnonCLC.txt
grep -Ff $Exinator1genes CISelements_nonCLC2_ENSGs.txt > overlapExinator1nonCLC.txt
grep -Ff $Exinator2genes CISelements_nonCLC2_ENSGs.txt > overlapExinator2nonCLC.txt

grep -Ff $PCAWGgenes CISelements_CLC2_ENSGs.txt > overlapPCAWGCLC.txt
grep -Ff $Exinator1genes CISelements_CLC2_ENSGs.txt > overlapExinator1CLC.txt
grep -Ff $Exinator2genes CISelements_CLC2_ENSGs.txt > overlapExinator2CLC.txt

echo "PCAWG/Exinator in hCIS hits"
cat overlapPCAWGnonCLC.txt 
cat overlapExinator1nonCLC.txt 
cat overlapExinator2nonCLC.txt 
cat overlapPCAWGCLC.txt 
cat overlapExinator1CLC.txt 
cat overlapExinator2CLC.txt 




#remove large files to save space

rm GENCODE.bed
rm GENCODE.txt
rm GENCODE1.txt
rm GENCODEall.bed
#rm GENCODEallids.bed
rm GENCODEallids.sorted.bed
rm GENCODEids.bed
rm GENCODEm1.bed
rm GENCODEnew.txt
rm proteinCodingGenes.txt


# make table for figure 6D with proportions


echo "Gene_type	Feature	Value" > OUTPUTforBarplot.txt

awk '{SUM += $3-$2} END {print "eintergenic\tAll_nucleotides\t"SUM}' correct_wholeGenome3.bed >> OUTPUTforBarplot.txt
awk '{SUM += $3-$2} END {print "cCGC\tAll_nucleotides\t"SUM}' CGCgenes1_M.bed >> OUTPUTforBarplot.txt
awk '{SUM += $3-$2} END {print "dnonCGC\tAll_nucleotides\t"SUM}' correct_nonCGC.bed >> OUTPUTforBarplot.txt
awk '{SUM += $3-$2} END {print "aCLC\tAll_nucleotides\t"SUM}' correct_CLC1.bed >> OUTPUTforBarplot.txt
awk '{SUM += $3-$2} END {print "bnonCLC\tAll_nucleotides\t"SUM}' correct_nonCLC2.bed >> OUTPUTforBarplot.txt

              
awk '{print $4}' CISelements_intergenic.bed | sort -u | wc -l | awk '{print "eintergenic\thCISelements\t"$1}' >> OUTPUTforBarplot.txt
awk '{print $4}' CISelements_CGC.bed | sort -u | wc -l | awk '{print "cCGC\thCISelements\t"$1}' >> OUTPUTforBarplot.txt
awk '{print $4}' CISelements_nonCGC.bed | sort -u | wc -l | awk '{print "dnonCGC\thCISelements\t"$1}' >> OUTPUTforBarplot.txt
awk '{print $4}' CISelements_CLC2.bed | sort -u | wc -l | awk '{print "aCLC\thCISelements\t"$1}' >> OUTPUTforBarplot.txt
awk '{print $4}' CISelements_nonCLC2.bed | sort -u | wc -l | awk '{print "bnonCLC\thCISelements\t"$1}' >> OUTPUTforBarplot.txt



Rscript /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/scripts/barplotFigure6.R

cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/plot1.pdf /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/plot1.pdf

cp /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/plot2.pdf /Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/$output_folder$2$3$4$5/results/plot2.pdf



cat $CLC2GENES | sort | uniq -c | awk '{print $2}' > CLCgenes_unique.txt
#true CLC genes 
grep -Ff CISelements_CLC2_ENSGs_lncRNAannotation.txt CLCgenes_unique.txt > CLCgenesWithCIS.txt


#genes that are not in CLC but overlay a CLC gene and share a CIS element
grep -vf CLCgenesWithCIS.txt CISelements_CLC2_ENSGs_lncRNAannotation.txt > lncRNAgenesOverlappingCLCgeneAndCIS.txt







#################################################################################
#################################################################################



# altered end November 2018

####new section to detect what region the CIS overlay that are in intergenic bin

#anaylze intergenic CIS in more detail, check which is closest gene to the CIS
#with -t report all genes, on both sides
#sort correctly!! 



sort -k1,1 -k2,2n CISelements_intergenic.bed | gsed 's/ /\t/g'  > CISelements_intergenicSorted.bed 
sort -k1,1 -k2,2n GENCODEallids.bed | cut -f1 -d"." | sed 's/"//g' | gsed 's/ /\t/g' > GENCODEallidsSorted.bed 


closestBed -d -a CISelements_intergenicSorted.bed -b GENCODEallidsSorted.bed > ClosestGeneIntergenic.bed
#closestBed -t first -a CISelements_intergenicSorted.bed -b GENCODEallidsSorted.bed > ClosestGeneIntergenic.bed


awk '{print $4}' CISelements_intergenic.bed | sort -u | wc -l
#209


awk '{print $4}' ClosestGeneIntergenic.bed | sort -u | wc -l
#209
awk '{print $4}' ClosestGeneIntergenic.bed | wc -l

awk '{print $5,$6,$7,$4,$8,$9}' ClosestGeneIntergenic.bed | gsed 's/ /\t/g' > ClosestGeneIntergenicSelection.bed
sort -k1,1 -k2,2n ClosestGeneIntergenicSelection.bed  > ClosestGeneIntergenicSorted.bed


#find in which categories the closest genes lies

intersectBed -a ClosestGeneIntergenicSorted.bed -b CGCgenes1_M.bed > closestGeneToCISelements_CGC.bed
intersectBed -v -a ClosestGeneIntergenicSorted.bed -b closestGeneToCISelements_CGC.bed > closestGeneToCISelements_afterCGC.bed 
#intersectBed -v -a outputCISVAR.bed -b CISelements_CGC.bed > CISelements_afterCGC.bed 


intersectBed -a closestGeneToCISelements_afterCGC.bed -b correct_nonCGC.bed > closestGeneToCISelements_nonCGC.bed
intersectBed -v -a closestGeneToCISelements_afterCGC.bed -b closestGeneToCISelements_nonCGC.bed > closestGeneToCISelements_afterNonCGC.bed 


intersectBed -a closestGeneToCISelements_afterNonCGC.bed -b correct_CLC1.bed > closestGeneToCISelements_CLC2.bed
intersectBed -v -a closestGeneToCISelements_afterNonCGC.bed -b closestGeneToCISelements_CLC2.bed > closestGeneToCISelements_afterCLC2.bed 


intersectBed -a closestGeneToCISelements_afterCLC2.bed -b correct_nonCLC2.bed > closestGeneToCISelements_nonCLC2.bed
intersectBed -v -a closestGeneToCISelements_afterCLC2.bed -b closestGeneToCISelements_nonCLC2.bed > closestGeneToCISelements_afterNonCLC2.bed 


intersectBed -a closestGeneToCISelements_afterNonCLC2.bed -b correct_wholeGenome3.bed > closestGeneToCISelements_intergenic.bed
intersectBed -v -a closestGeneToCISelements_afterNonCLC2.bed -b closestGeneToCISelements_intergenic.bed > closestGeneToCISelements_afterIntergenic.bed 
cat closestGeneToCISelements_afterIntergenic.bed | wc -l



#how many of the genes are already in the nonCLC pot?
awk '{print $5}' closestGeneToCISelements_nonCLC2.bed | sort | uniq -c | awk '{print $2}' >  closestGeneToCISelements_nonCLC2_ENSGs.txt


grep -Ff closestGeneToCISelements_nonCLC2_ENSGs.txt CISelements_nonCLC2_ENSGs_lncRNA.bed | awk '{print $8}' | sort | uniq -c | awk '{print $2}'  > nonCLCENSGsInLiftOverAndClosests.txt

grep -vf  nonCLCENSGsInLiftOverAndClosests.txt closestGeneToCISelements_nonCLC2_ENSGs.txt  > nonCLCENSGsOnlyInClosests.txt

#16 together, correct!







############

############

############

############ new section, overlay intergenic CIS with miTranscriptome to see if they overlay with non GENCODE id transcripts


sort -k1,1 -k2,2n CISelements_intergenic.bed | gsed 's/ /\t/g'  > CISelements_intergenicSorted.bed 

miTranscriptome="/Users/ninavan/Dropbox/CLC2analysis/CISanalysis/LiftOverResults/MiTranscriptome/mitranscriptome.bed"


awk '{print $1,$2,$3,$4}' $miTranscriptome | sort | gsed 's/ /\t/g'  >  miTranscriptomeGenomicLocations.bed


intersectBed -wa -wb -a CISelements_intergenicSorted.bed -b miTranscriptomeGenomicLocations.bed > CISelements_intergenic_miTranscriptome.bed

#check how many intergenic CIS have a transcript in miTranscriptome
awk '{print $4}' CISelements_intergenic_miTranscriptome.bed | sort | uniq -c | wc -l
#186 from 229 in total CISelements_intergenicSorted.bed


############
############
############
############
############
############


#many genes intersect with transcripts so maybe use different database, one more stringent than miTranscriptome









#further analyze the CIS in CLC and nonCLC, and show in what cancer type they are present



cat CISelementsALLINFO_nonCLC_FINAL_lncRNAannotation.txt | sort | uniq -c | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14}' > CISelementsALLINFO_nonCLC_lncRNAannotation_Unique.txt 




#cat CISelementsALLINFO_nonCLC_lncRNAannotation.txt | sort | uniq -c | awk '{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' > CISelementsALLINFO_nonCLC_lncRNAannotation_Unique.txt 

awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' CISelementsALLINFO_nonCLC_lncRNAannotation_Unique.txt | sort  >  CISelementsALLINFO_nonCLC_lncRNAannotation_Unique_SortedBasedOnENSG.txt

awk '$12 == "BloodCancer"' CISelementsALLINFO_nonCLC_lncRNAannotation_Unique_SortedBasedOnENSG.txt > nonCLC_CISelements_inBloodCancer.txt

#check this to see if all are unique, in this case yes!
#but there are CIS that map to two loci 
awk '{print $10,$1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$12,$13}' nonCLC_CISelements_inBloodCancer.txt | sort | uniq -c 




#find in all cancers 


awk '{print $1, $2, $12, $13}' CISelementsALLINFO_nonCLC_lncRNAannotation_Unique_SortedBasedOnENSG.txt


# | sort | uniq -c > nonCLC_CISelements_inBloodCancer.txt

#check this to see if all are unique, in this case yes!
#but there are CIS that map to two loci 
awk '{print $10,$1,$2,$3,$4,$5,$6,$7,$8,$9,$11,$12,$13}' nonCLC_CISelements_inBloodCancer.txt | sort | uniq -c 




awk '{print $1}' CISelementsALLINFO_nonCLC_lncRNAannotation_Unique_SortedBasedOnENSG.txt | wc -l





####################################################################################

echo "end of script"

####################################################################################
#end
####################################################################################



#incorporate conservation score to each gene and espression data from UCSC






# to do list

# also sort out intergenic hCIS hits that map to mouse protein coding

# redo figures with correct input! what to do with the ones sorted as mouse protein coding?




























