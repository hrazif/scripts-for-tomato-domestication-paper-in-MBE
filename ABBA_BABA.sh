#!/bin/bash

# nohup ./job6.sh >>job6.output &



 for i in 01 02 03 04 05 06 07 08 09 10 11 12; do

#input split into separate chromosomes
input=ch"$i".geno.gz

#-f phased; it doesn't really mean phased; -w: window size in bp
 python /home/caicedo/Hamid/20180813_redoing_analyses_based_on_new_tree/6_ABBA_BABA_using_Simons_new_scripts/test_with_Simons_scripts/genomics_general/freq.py -g $input \
 -p SP_SECU -p SP_PER -p SP_NECU -p SLC_ECU -p SLC_PER \
 -p SLC_San_Martin -p SLC_MEX_CA_NSA -p SLC_MEX -p SLL -p SPENN \
 --popsFile populations.txt --target derived \
 -o chr"$i".derFreq.tsv.gz &


 done


 R CMD BATCH ABBA_BABA.r
