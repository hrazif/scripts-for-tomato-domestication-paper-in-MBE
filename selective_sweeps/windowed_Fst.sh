#!/bin/bash 

#nohup ./job10.sh>>job10.output&

#input identifical to fastStructure input, except not filtering SNPs based on missing data
VCF_Xinvariant=/home/caicedo/Hamid/20170612/1_selctive_sweeps_all_SP_and_SLC_ECU/InitFilt_exWilds_exF1_xInd25p_xProblematic_Xinvariant.vcf.gz

#for domestication
~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --weir-fst-pop SP.txt --weir-fst-pop SLC_ECU.txt --fst-window-size 10000  --maf 0.025 --out 20190311_3_10_SP_to_SLC_ECU&
~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --weir-fst-pop SP.txt --weir-fst-pop SLC_ECU_PER_SM.txt --fst-window-size 10000  --maf 0.025 --out 20190311_3_10_SP_to_SLC_ECU_PER_SM&

~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --weir-fst-pop SP_NECU.txt --weir-fst-pop SLC_ECU.txt --fst-window-size 10000  --maf 0.025 --out 20190311_3_10_SP_N_ECU_to_SLC_ECU&
~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --weir-fst-pop SP_NECU.txt --weir-fst-pop SLC_ECU_PER_SM.txt --fst-window-size 10000  --maf 0.025 --out 20190311_3_10_SP_N_ECU_to_SLC_ECU_PER_SM&


~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --weir-fst-pop SLC_MEX.txt --weir-fst-pop SLL.txt --fst-window-size 10000  --maf 0.025 --out 20190311_3_10_SLC_MEX_to_SLL&

#for dedomestication
~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --weir-fst-pop SLC_ECU_PER_SM.txt --weir-fst-pop SLC_MEX_CA_NSA.txt --fst-window-size 10000  --maf 0.025 --out 20190311_3_10_SLC_ECU_PER_SM_to_SLC_MEX_CA_NSA&
~/usr/bin/vcftools --gzvcf $VCF_Xinvariant --weir-fst-pop SLC_ECU_PER_SM.txt --weir-fst-pop SLC_MEX.txt --fst-window-size 10000  --maf 0.025 --out 20190311_3_10_SLC_ECU_PER_SM_to_SLC_MEX&

