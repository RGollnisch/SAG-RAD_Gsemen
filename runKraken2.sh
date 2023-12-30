# Run kraken2 from the cleaned reads
# Specify the directory of the processed fastq files
FQDIR=.

#Specify the number of threads to be used
THREADS=100

for i in $FQDIR/*gz; do kraken2 -db download-taxonomy --threads $THREADS --use-names  --gzip-compressed -output ${i#}.krakenout $i; echo "${i#}.krakenout";  done

