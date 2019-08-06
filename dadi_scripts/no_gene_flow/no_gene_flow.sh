#!/bin/bash


#nohup ./job7_1.sh >> job7_1.out &

#input identical with the fastStructure input, except including S. penelli (outgroup) to polarize the SNPs.
zcat /home/caicedo/Hamid/20180813_redoing_analyses_based_on_new_tree/2_dadi_dating/8_full_genome_SP_NECU_and_SLCSLL/no_gene_flow/structure_input_wOutgroup_simple_xChr0.vcf.gz >  structure_input_wOutgroup_simple_xChr0.vcf


mv structure_input_wOutgroup_simple_xChr0.vcf full_genome.vcf

perl ~/software/convert_vcf_to_dadi_input.pl full_genome.vcf pops.txt

rm full_genome.vcf
mv full_genome.vcf.data full_genome.dadi



mkdir 2d_comp_fig 2d_sfs_fig out
python run_dadi_f2pops_split_no_mig.py