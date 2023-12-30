#!/bin/bash

# Run a loop to run one replicate for each K value
# Samples from USA has 10 populations, hence K 2-10 is run
# By changing the rand value in the -o parameter (rand1, rand2 ..) the output files will not overwrite previous replication.

# Check the number of variant sites retained (last value) and modify the -L value accordingly. 
#grep Kept populations.log

#Check the number of samples & populations and adjust -L for samples and K 2 to populations.
#grep Working populations.log

# variable r= Number of replications and variable k= the population sizes to be tested. Note that r*k jobs will be run in parallel.

#Go to structure output directory and remove old structure files 
cd ./structure
rm populations.structure
rm structure.EDT.txt

#Specify the file for extraparams in Structure
#extraparams=extraparams
extraparams=extraparams_FREQSCORR_0

#Modify the Stack structure output to fit Structure

# USA samples
#
#populations_p=7
#populations_r=0.7

# USA minus SM & LL populations 
populations_p=7
populations_r=0.7
cp ../populations/popmap_50k_usa_minus_SM_LL/p${populations_p}r${populations_r}/populations.structure .
echo "Copying Populations results from ../populations/popmap_50k_usa_minus_SM_LL/p${populations_p}r${populations_r}/populations.structure  as Structure input"

# Generate a Structure file based on the naming of the populations 
# Needs to be changed to fit the populations of the current study

cat populations.structure| grep -v "^# Stacks"|sed 's/SM	/1	/'| sed 's/CC	/2	/'|sed 's/CS	/3	/'|sed 's/LL	/4	/'|sed 's/NO	/5	/'|sed 's/PQ	/6	/'|sed 's/PT	/7	/'|sed 's/RE	/8	/'|sed 's/VO	/9	/'|sed 's/WN	/10	/' | sed 's/BR	/11	/'| sed 's/GF	/12	/'|sed 's/HE	/13	/'|sed 's/HU	/14	/'|sed 's/KO	/15	/'|sed 's/NT	/16	/'|sed 's/PB	/17	/'|sed 's/PE	/18	/'|sed 's/PI	/19	/'|sed 's/PS	/20	/'|sed 's/RM	/21	/'|sed 's/RO	/22	/'|sed 's/SE	/23	/'|sed 's/SI	/24	/'|sed 's/VP	/25	/'> structure.EDT.txt

#Settings for USA samples minus SM & LL
RVALUE=30
KVALUE=8
Samples=138
Variants=715
Continent=USA
INPUTFILE=structure.EDT.txt

#Loop through the r and k values to get Structure results for each combination

for r in $(seq 21 ${RVALUE});
  do
    for k in $(seq 2 ${KVALUE});
      do
        # Make outdirs for every k tested.
        mkdir -p ./${Continent}/K${k}; 
        structure -e ${extraparams} -D $RANDOM -K ${k} -L ${Variants} -N ${Samples} -i $INPUTFILE -o ./${Continent}/K${k}/${Continent}_${populations_p}.${populations_r}.${r}_k${k}& 
      done;
  done


