#!/bin/bash 

#nohup ./job2_all_pops.sh>>job2_all_pops.output&

#input file same as FastStructure input, except without filtering for missing SNPs
input_vcf=/home/caicedo/Hamid/20170612/1_selctive_sweeps_all_SP_and_SLC_ECU/InitFilt_exWilds_exF1_xInd25p_xProblematic_Xinvariant.vcf.gz

for i in SP_SECU SP_PER SP_NECU SLC_ECU SLC_PER SLC_San_Martin SLC_MEX_CA_NSA SLC_MEX SLL; do


	mkdir "$i"
	cd "$i"
	cp ../"$i"_plink.txt ./
	for j in 01 02 03 04 05 06 07 08 09 10 11 12; do
		#--ld-window 1000000 is the same as the length of the LD window (1000 kb); keeping all SNPs; --bp-space (thin in vcftools) to keep one snp per 1 kb
		~/usr/bin/plink --allow-extra-chr --vcf-half-call m  --vcf $input_vcf --keep "$i"_plink.txt --maf 0.001 --bp-space 1000 --chr SL2.50ch"$j" --r2 gz --ld-window-kb 1000 --ld-window-r2 0 --ld-window 1000000 --out "$i".chr"$j" --threads 20
		gzip -cd "$i".chr"$j".ld.gz|grep -v "BP"|nawk -F " " '{print ($5-$2)/1000 " " $7}'|awk '{printf "%.0f %.5f\n", $1, $2}'| gzip -c > "$i".chr"$j".simple.ld.gz
		rm "$i".chr"$j".ld.gz
	done


	for k in $(seq 1 1000); do zcat "$i".chr*.simple.ld.gz|awk '$1 == '$k''|datamash -W -g 1 mean 2 >> distance_meanR2;done &

	cd .. 
done



