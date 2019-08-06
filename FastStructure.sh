#!/bin/bash 

#bsub -q sgi -W 240:00 -n 1 -R rusage[mem=40000] -o job11.out -e job11.err ./job11.sh

module load python3/3.5.0 fastStructure/v1.0-4-ge47212f R/3.3.1

vcftools --gzvcf /project/uma_ana_caicedo/Hamid/20170220/1/newVCF_withIntialFiltering.vcf.gz --remove wilds_F1.txt --recode -c | gzip -c > InitFilt_exWilds_exF1.vcf.gz

vcftools --gzvcf /project/uma_ana_caicedo/Hamid/20170410/1_fastSRUCTURE_combined_dataset_indv25p_SNPs10p/InitFilt_exWilds_exF1_xInd25p.vcf.gz --remove problematic.txt --recode -c | gzip -c > InitFilt_exWilds_exF1_xInd25p_xProblematic.vcf.gz

#--maf for removing singletons; --max-missing-count for removing SNPs missing in > 10% of individuals
vcftools --gzvcf InitFilt_exWilds_exF1_xInd25p_xProblematic.vcf.gz --max-alleles 2 --maf 0.00680 --max-missing-count 30 --recode -c | gzip -c > InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_biallelic_xSingles.vcf.gz

vcftools --gzvcf InitFilt_exWilds_exF1_xInd25p_xProblematic.vcf.gz --min-alleles 3 --max-missing-count 30 --recode -c | gzip -c > InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_multiallelic.vcf.gz

#to extract the list of biallelic and multi-allelic SNPs and combine them
zcat InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_biallelic_xSingles.vcf.gz|grep -v "#"|cut -f3>InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_biallelic_xSingles_SNPs
zcat InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_multiallelic.vcf.gz|grep -v "#"|cut -f3>InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_multiallelic_SNPs
cat InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_biallelic_xSingles_SNPs InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_multiallelic_SNPs|sort -V> InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_all_SNPs

#to keep only multi-allelic and biallelic SNPs (that were filtered for singles)
vcftools --gzvcf InitFilt_exWilds_exF1_xInd25p_xProblematic.vcf.gz --snps InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_all_SNPs --recode -c | gzip -c > InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_all.vcf.gz

#to make plink files to use as input for fastStructure
/home/hr42a/software/plink/1.09/plink  --allow-extra-chr --vcf-half-call m  --vcf  InitFilt_exWilds_exF1_xInd25p_xProblematic_xSNPs10p_all.vcf.gz

#running FastStructure
for i in {2..10};
do 
bsub -q long -W 480:00  -R rusage[mem=80000]  -n 1 -o job11_$i.out -e job11_$i.err python /share/pkg/fastStructure/v1.0-4-ge47212f/structure.py -K $i --input=plink --output=20170417_11 --full; 

done

#chooseK algorithm to choose K that fits the data best; run this after fastStructure finishes
##python /share/pkg/fastStructure/v1.0-4-ge47212f/chooseK.py --input=20170417_11
