# Set the directory for sstacks and popmap 
cd ./sstacks/continents25kQ2_catalog200

PopMap=./ustacks/info/popmap_sstacks25k

sstacks -P ./ -M $PopMap -p 200
