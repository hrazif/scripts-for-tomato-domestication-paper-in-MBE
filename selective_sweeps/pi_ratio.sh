#!/bin/bash 

#nohup ./job10.sh>>job10.output&

#input identifical to fastStructure input, except not filtering SNPs based on missing data
VCF_Xinvariant=/home/caicedo/Hamid/20170612/1_selctive_sweeps_all_SP_and_SLC_ECU/InitFilt_exWilds_exF1_xInd25p_xProblematic_Xinvariant.vcf.gz

for i in SP SLC_ECU SLC_ECU_PER_SM SLC_MEX_CA_NSA SLC_MEX SLL; do
~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --keep "$i".txt  --window-pi 10000  --out "$i"&
done
