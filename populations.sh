#!/bin/bash

#Specify the sstacks directory
DIR=./sstacks

cd $DIR

#---------------Run both continents ------------------------
# Specify the directories for popmap, output and the p, r parameters in Stacks

PopMap="../4.ustacks/info/popmap50k_continents_minus_SM_LL_20pop"
PopOut=popmap_50k_continents_minus_SM_LL_20pop
p=20
r="0.6"
#-----------------------------------------------------------

mkdir -p  ../populations/${PopOut}/p${p}r${r}
populations -P ./ --popmap ${PopMap} -O ../populations/${PopOut}/p${p}r${r} -p ${p} -r ${r} -f p_value -t 60 --structure --genepop --write-single-snp --fstats 

#----------------------------------------------------------

