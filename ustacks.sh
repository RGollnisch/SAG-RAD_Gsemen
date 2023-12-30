#!/bin/bash -l
#SBATCH --job-name=ustacks_lib01-A
#SBATCH --account snic2017-7-349
#SBATCH -p core -n 16
#SBATCH -t 12:00:00

#module load bioinfo-tools Stacks/2.59

# Make sure that you are running the conda environment using
#  conda activate PopGen

	IN_DIRS=(
UE-2976-lib01-A
UE-2976-lib02-B
UE-2976-lib03-C
UE-2976-lib04-D
UE-2976-lib05-E
UE-2976-lib06-F
UE-2976-lib07-A
UE-2976-lib08-B
UE-2976-lib09-C
UE-2976-lib10-D
UE-2976-lib11-E
UE-2976-lib12-F
UE-2976-lib13-A
UE-2976-lib14-B
UE-2976-lib15-C
UE-2976-lib16-D
UE-2976-lib17-E
UE-2976-lib18-F
)

for dir in ${IN_DIRS[@]}
do
	SUBDIR=$dir

	INFO=($(echo $SUBDIR|sed 's/UE-2976-//'))

(${CLEAN_DIR##*/})
	M_VALUE=3
	LT_m_VALUE=5

	CURR_DIR=./cleaned/${SUBDIR}
	DENOVO_DIR=./ustacks/${M_VALUE}/${LT_m_VALUE}/${SUBDIR}
	ID_FILE=./ustacks/info/ustacks_IDs_${INFO}.csv

	# echo "$CURR_DIR is current dir, Parametenrs M=${M_VALUE} m=${LT_m_VALUE}"
	mkdir -p $DENOVO_DIR

	file=($(cut -d',' -f4 $ID_FILE))
	name=($(cut -d',' -f5 $ID_FILE))
	ustacks_id=($(cut -d',' -f6 $ID_FILE))

	file=(${file[@]:1})
	name=(${name[@]:1})
	ustacks_id=(${ustacks_id[@]:1})

	# default ustacks: M=2, m=3, N=M+2
	# gonyPop ustacks: M=3, m=5, N=1

# Loop to run ustacks for all fastq file. 
# Note that the "&" at the end of the ustacks command will make the commands run in parallel. Remove to reduce the computer requirements at the expense of lower throughput

	for i in  ${!file[@]}; do
        	ustacks -f ${CURR_DIR}/${file[$i]} -o $DENOVO_DIR -i ${ustacks_id[$i]} --name ${name[$i]} -M ${M_VALUE} -m ${LT_m_VALUE} -p 16 > ustacks.${name[$i]}_M${M_VALUE}_m${LT_m_VALUE}.final.out &
	done
done
