#!/bin/bash

#nohup ./job11.sh>>job11.output&

#input identifical to fastStructure input, except not filtering SNPs based on missing data
VCF=/home/caicedo/Hamid/20170612/1_selctive_sweeps_all_SP_and_SLC_ECU/InitFilt_exWilds_exF1_xInd25p_xProblematic_Xinvariant.vcf.gz

#--maf 0.001 for removing invariants
#--max-maf for keeeping rare alleles

#to create a list of rare alleles in SLL
~/usr/bin/vcftools --gzvcf $VCF --max-alleles 2 --keep SLL.txt --maf 0.001 --max-maf 0.05 --extract-FORMAT-info GT --out "SLL_rare_alleles"
#removing potentially admixed accessions
~/usr/bin/vcftools --gzvcf $VCF --max-alleles 2 --keep SLL_xAdmixed.txt --maf 0.001 --max-maf 0.05 --extract-FORMAT-info GT --out "SLL_xAdmixed_rare_alleles"
cut -f1-2 SLL_xAdmixed_rare_alleles.GT.FORMAT|sed '1d'> SLL_xAdmixed_rare_alleles_pos



#to extract remaining SLL rare alleles in SLC_MEX
~/usr/bin/vcftools --gzvcf $VCF --max-alleles 2 --maf 0.001 --keep SLC_MEX.txt --positions SLL_xAdmixed_rare_alleles_pos --extract-FORMAT-info GT  --out "SLL_xAdmixed_rare_alleles_in_SLC_MEX"
cut -f1-2 SLL_xAdmixed_rare_alleles_in_SLC_MEX.GT.FORMAT|sed '1d'> SLL_xAdmixed_rare_alleles_in_SLC_MEX_pos
#to make the list of positions for the rare SLL_xAdmixed alleles not present in SLC_MEX
grep -vwF -f SLL_xAdmixed_rare_alleles_in_SLC_MEX_pos SLL_xAdmixed_rare_alleles_pos >SLL_xAdmixed_rare_alleles_not_in_SLC_MEX_pos


#to extract remaining SLL rare alleles in SLC_MEX_CA_NSA
~/usr/bin/vcftools --gzvcf $VCF --max-alleles 2 --maf 0.001 --keep SLC_MEX_CA_NSA.txt --positions SLL_xAdmixed_rare_alleles_not_in_SLC_MEX_pos --extract-FORMAT-info GT  --out "SLL_xAdmixed_rare_alleles_in_SLC_MEX_CA_NSA"
cut -f1-2 SLL_xAdmixed_rare_alleles_in_SLC_MEX_CA_NSA.GT.FORMAT|sed '1d'> SLL_xAdmixed_rare_alleles_in_SLC_MEX_CA_NSA_pos
#to make the list of positions for the rare SLL_xAdmixed alleles not present in SLC_MEX_CA_NSA
grep -vwF -f SLL_xAdmixed_rare_alleles_in_SLC_MEX_CA_NSA_pos SLL_xAdmixed_rare_alleles_not_in_SLC_MEX_pos >SLL_xAdmixed_rare_alleles_not_in_SLC_MEX_CA_NSA_pos


#to extract remaining SLL rare alleles in SLC_San_Martin
~/usr/bin/vcftools --gzvcf $VCF --max-alleles 2 --maf 0.001 --keep SLC_San_Martin.txt --positions SLL_xAdmixed_rare_alleles_not_in_SLC_MEX_CA_NSA_pos --extract-FORMAT-info GT  --out "SLL_xAdmixed_rare_alleles_in_SLC_San_Martin"
cut -f1-2 SLL_xAdmixed_rare_alleles_in_SLC_San_Martin.GT.FORMAT|sed '1d'> SLL_xAdmixed_rare_alleles_in_SLC_San_Martin_pos
#to make the list of positions for the rare SLL_xAdmixed alleles not present in SLC_San_Martin
grep -vwF -f SLL_xAdmixed_rare_alleles_in_SLC_San_Martin_pos SLL_xAdmixed_rare_alleles_not_in_SLC_MEX_CA_NSA_pos >SLL_xAdmixed_rare_alleles_not_in_SLC_San_Martin_pos

#to extract remaining SLL rare alleles in SLC_PER
~/usr/bin/vcftools --gzvcf $VCF --max-alleles 2 --maf 0.001 --keep SLC_PER.txt --positions SLL_xAdmixed_rare_alleles_not_in_San_Martin_pos --extract-FORMAT-info GT  --out "SLL_xAdmixed_rare_alleles_in_SLC_PER"
cut -f1-2 SLL_xAdmixed_rare_alleles_in_SLC_PER.GT.FORMAT|sed '1d'> SLL_xAdmixed_rare_alleles_in_SLC_PER_pos
#to make the list of positions for the rare SLL_xAdmixed alleles not present in SLC_PER
grep -vwF -f SLL_xAdmixed_rare_alleles_in_SLC_PER_pos SLL_xAdmixed_rare_alleles_not_in_San_Martin_pos >SLL_xAdmixed_rare_alleles_not_in_SLC_PER_pos

#to extract remaining SLL rare alleles in SLC_ECU
~/usr/bin/vcftools --gzvcf $VCF --max-alleles 2 --maf 0.001 --keep SLC_ECU.txt --positions SLL_xAdmixed_rare_alleles_not_in_SLC_PER_pos --extract-FORMAT-info GT  --out "SLL_xAdmixed_rare_alleles_in_SLC_ECU"
cut -f1-2 SLL_xAdmixed_rare_alleles_in_SLC_ECU.GT.FORMAT|sed '1d'> SLL_xAdmixed_rare_alleles_in_SLC_ECU_pos
#to make the list of positions for the rare SLL_xAdmixed alleles not present in SLC_ECU
grep -vwF -f SLL_xAdmixed_rare_alleles_in_SLC_ECU_pos SLL_xAdmixed_rare_alleles_not_in_SLC_PER_pos >SLL_xAdmixed_rare_alleles_not_in_SLC_ECU_pos

#to extract remaining SLL rare alleles in SP_NECU
~/usr/bin/vcftools --gzvcf $VCF --max-alleles 2 --maf 0.001 --keep SP_NECU.txt --positions SLL_xAdmixed_rare_alleles_not_in_SLC_ECU_pos --extract-FORMAT-info GT  --out "SLL_xAdmixed_rare_alleles_in_SP_NECU"
cut -f1-2 SLL_xAdmixed_rare_alleles_in_SP_NECU.GT.FORMAT|sed '1d'> SLL_xAdmixed_rare_alleles_in_SP_NECU_pos
#to make the list of positions for the rare SLL_xAdmixed alleles not present in SP_NECU
grep -vwF -f SLL_xAdmixed_rare_alleles_in_SP_NECU_pos SLL_xAdmixed_rare_alleles_not_in_SLC_ECU_pos >SLL_xAdmixed_rare_alleles_not_in_SP_NECU_pos

#to extract remaining SLL rare alleles in SP_PER
~/usr/bin/vcftools --gzvcf $VCF --max-alleles 2 --maf 0.001 --keep SP_PER.txt --positions SLL_xAdmixed_rare_alleles_not_in_SP_NECU_pos --extract-FORMAT-info GT  --out "SLL_xAdmixed_rare_alleles_in_SP_PER"
cut -f1-2 SLL_xAdmixed_rare_alleles_in_SP_PER.GT.FORMAT|sed '1d'> SLL_xAdmixed_rare_alleles_in_SP_PER_pos
#to make the list of positions for the rare SLL_xAdmixed alleles not present in SP_PER
grep -vwF -f SLL_xAdmixed_rare_alleles_in_SP_PER_pos SLL_xAdmixed_rare_alleles_not_in_SP_NECU_pos >SLL_xAdmixed_rare_alleles_not_in_SP_PER_pos

#to extract remaining SLL rare alleles in SP_SECU
~/usr/bin/vcftools --gzvcf $VCF --max-alleles 2 --maf 0.001 --keep SP_SECU.txt --positions SLL_xAdmixed_rare_alleles_not_in_SP_PER_pos --extract-FORMAT-info GT  --out "SLL_xAdmixed_rare_alleles_in_SP_SECU"
cut -f1-2 SLL_xAdmixed_rare_alleles_in_SP_SECU.GT.FORMAT|sed '1d'> SLL_xAdmixed_rare_alleles_in_SP_SECU_pos
#to make the list of positions for the rare SLL_xAdmixed alleles not present in SP_SECU
grep -vwF -f SLL_xAdmixed_rare_alleles_in_SP_SECU_pos SLL_xAdmixed_rare_alleles_not_in_SP_PER_pos >SLL_xAdmixed_rare_alleles_not_in_SP_SECU_pos



