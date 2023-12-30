
#Set the directory of sstacks
DIR=./sstacks/continents25kQ2_catalog200/
cd $DIR

PopMap="../../ustacks/info/popmap_sstacks25k"
ReadsDir="../../kraken2/Unclassified/"

#ReadsDir should point to a directory with the filtered (unclassified according to Kraken2) paired end reads are located.
#ReadsDir="../2.kraken2/Unclassified/"

tsv2bam -P ./ -M ./${PopMap} -t 100 -R ${ReadsDir}
