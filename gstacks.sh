#Specify the sstacks and popmap catalog
DIR=./sstacks/continents25kQ2_catalog200/
PopMap="popmap_sstacks25k"

#Specify directory for the Kraken processed fastq files
ReadsDir="./kraken2/Unclassified/"

cd $DIR
gstacks -P ./ -M ./${PopMap} -t 100 --rm-pcr-duplicates
