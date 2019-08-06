#!/bin/bash 

#nohup ./job3.sh>>job3.output &


#to convert the fastStructure input to TreeMix input, excluding unassigned accessions from data.clust
~/usr/bin/plink --allow-extra-chr --bfile '/home/caicedo/Hamid/20170417/11_fastSTRUCTURE_like_6_removing_more_problematic_accessions [good]'/plink --freq --within data.clust 

gzip plink.frq.strat
python ~/software/plink2treemix.py plink.frq.strat.gz treemix.frq.gz


#running TreeMix; -k 12000, representing 10 MB as in the TreeMix paper
for i in {0..10}
do ~/usr/bin/treemix -i treemix.frq.gz -root SP_SECU -k 12000 -m $i -se -o 20190311_3."$i" &
done


##run the codes below after TreeMix finishes running to summarize the results.
##grep 'Exiting' *.llik|sed s/'20190311_3.'//g|sed s/".llik:Exiting ln(likelihood) with [0-9]* migration events: "/"\t"/g > likelihoods.txt

#for plotting the TreeMix results; plotting_funcs.R is provided in the TreeMix package
##R CMD BATCH plotting_funcs.R