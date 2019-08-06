#!/bin/bash

# nohup ./job2_2.sh > job2_2.out & 

#input identical with the fastStructure input, except including S. penelli (outgroup) to polarize the SNPs.
cp /home/caicedo/Hamid/20190311_redoing_analyses_for_the_paper/7_dadi_dating/1_full_genome_SP_NECU_and_SLCSLL/no_gene_flow/full_genome.dadi ./

mkdir 2d_comp_fig 2d_sfs_fig out
python run_dadi_f2pops_split_with_mig.py