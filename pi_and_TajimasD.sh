#!/bin/bash 

#nohup ./job5.sh>>job5.output&

#input identifical to fastStructure input, except not filtering SNPs based on missing data
VCF_Xinvariant=/home/caicedo/Hamid/20170612/1_selctive_sweeps_all_SP_and_SLC_ECU/InitFilt_exWilds_exF1_xInd25p_xProblematic_Xinvariant.vcf.gz

for i in SP_SECU SP_PER SP_NECU SLC_ECU SLC_PER SLC_San_Martin SLC_MEX_CA_NSA SLC_MEX SLL SLL_xAdmixed; do
~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --keep "$i".txt --maf 0.001  --window-pi 10000 --out "$i" &

~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --keep "$i".txt --maf 0.001 --TajimaD 10000 --out "$i" &

done
