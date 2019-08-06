#!/bin/bash 

#nohup ./job10.sh>>job10.output&

input_vcf=/home/caicedo/Hamid/20190311_redoing_analyses_for_the_paper/2_LD/bigVCF_after_qual_filt.vcf.gz


#creating inputs for each group and including S. pennellii (outgroup) for polarizing the SNPs
 for i in SP SLC_ECU SLC_ECU_PER_SM SLC_MEX_CA_NSA SLC_MEX SLL; do
 ~/usr/bin/vcftools --gzvcf $input_vcf --keep "$i"_wOutgroup.txt --maf 0.001  --recode -c |sed 's/\<SL2.50ch01\>/1/g; s/\<SL2.50ch02\>/2/g; s/\<SL2.50ch03\>/3/g; s/\<SL2.50ch04\>/4/g; s/\<SL2.50ch05\>/5/g; s/\<SL2.50ch06\>/6/g; s/\<SL2.50ch07\>/7/g; s/\<SL2.50ch08\>/8/g; s/\<SL2.50ch09\>/9/g; s/\<SL2.50ch10\>/10/g; s/\<SL2.50ch11\>/11/g; s/\<SL2.50ch12\>/12/g'|bgzip -c > "$i".vcf.gz
 done



#splitting into separte chromosomes
for i in SP SLC_ECU SLC_ECU_PER_SM SLC_MEX_CA_NSA SLC_MEX SLL; do
	for j in {1..12}; do
		~/usr/bin/vcftools --gzvcf "$i".vcf.gz  --chr "$j"  --recode -c |bgzip -c > "$i"_chr"$j".vcf.gz 
 	done
done



#making SweeD input

for i in SP SLC_ECU SLC_ECU_PER_SM SLC_MEX_CA_NSA SLC_MEX SLL; do
		for j in {1..12}; do
			zcat "$i"_chr"$j".vcf.gz > "$i"_chr"$j".vcf
			perl ~/software/convert_vcf_to_dadi_input.pl "$i"_chr"$j".vcf "$i"_dadi_pop.txt
		    mv "$i"_chr"$j".vcf.data "$i"_chr"$j".dadi
		done	
done		
	

#running SweeD

for i in SP SLC_ECU SLC_ECU_PER_SM SLC_MEX_CA_NSA SLC_MEX SLL; do
	for j in {1..12}; do
		zcat "$i"_chr"$j".vcf.gz > "$i"_chr"$j".vcf	
		~/software/sweed/SweeD-P -name "$i"_chr"$j" -input "$i"_chr"$j".sweeD -grid 100000  -strictPolymorphic -noSeparator -threads 3
	 done &
done






