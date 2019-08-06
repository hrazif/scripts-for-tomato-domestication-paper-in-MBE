#!/bin/bash 

#./job5.sh

#input identifical to fastStructure input, except not filtering SNPs based on missing data
VCF_Xinvariant=/home/caicedo/Hamid/20170612/1_selctive_sweeps_all_SP_and_SLC_ECU/InitFilt_exWilds_exF1_xInd25p_xProblematic_Xinvariant.vcf.gz

for i in SP_SECU SP_PER SP_NECU SLC_ECU SLC_PER SLC_San_Martin SLC_MEX_CA_NSA SLC_MEX SLL SLL_xAdmixed; do
nohup ~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --keep "$i".txt --maf 0.001  --site-pi --out "$i" > "$i"_site_pi.out &
nohup ~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --keep "$i".txt --maf 0.001  --hardy --out "$i" > "$i"_hwe.out &
done

