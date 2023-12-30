# cstacks -P path_to_ustacksDIR -M popmap -n number_of_mismatches_allowed_between_sample_loci_in_catalog -p threads -o outDIR -c pregenerated_catalog

# Specify the directory where ustacks files are located
ustacksDIR="./ustacks/3/5"
# If popmap file is called popmap_cstacks25 then use popmap_cstacks and omit the 25, 50, 75 and 100.
CstacksDIR="./cstacks"
popMap="popmap_cstacksQ2_"
n_PARAM=3


cd ${CstacksDIR}

# Generate directories for the incrementally larger catalogs
mkdir -p ./ ./catalogQ2_25/ ./catalogQ2_50/ ./catalogQ2_75/ ./catalogQ2_100/ ./catalogQ2_125 ./catalogQ2_150 ./catalogQ2_175 ./catalogQ2_200

# --------------- First run without any catalog -------------------------------------

# Took 139 min to run first 25
time cstacks -P . -M ${popMap}25 -n ${n_PARAM} -p 200 >./catalogQ2_25/cstacks25.log
# Move catalog to catalog dir.
mv catalog*gz ./catalogQ2_25/

#----------------- Rerun cstacks with another 25 samples ----------------------------
#Second round using existing catalog
time cstacks -P . -M ${popMap}50 -n ${n_PARAM} -p 200 --catalog ./catalogQ2_25/catalog 
# Move catalog to catalog dir.
mv catalog*gz ./catalogQ2_50/

#----------------- Rerun cstacks with another 25 samples ----------------------------
# Third round
time cstacks -P . -M ${popMap}75 -n ${n_PARAM} -p 200 --catalog ./catalogQ2_50/catalog 
# Move catalog to catalog dir.
mv catalog*gz ./catalogQ2_75/

#----------------- Rerun cstacks with another 25 samples ----------------------------
#Fourth round 
time cstacks -P . -M ${popMap}100 -n ${n_PARAM} -p 200  --catalog ./catalogQ2_75/catalog 
# Move catalog to catalog dir.
mv catalog*gz ./catalogQ2_100/

#----------------- Rerun cstacks with another 25 samples ----------------------------
#Fifth round
time cstacks -P . -M ${popMap}125 -n ${n_PARAM} -p 200  --catalog ./catalogQ2_100/catalog
# Move catalog to catalog dir.
mv catalog*gz ./catalogQ2_125/


#----------------- Rerun cstacks with another 25 samples ----------------------------
#Sixth round
time cstacks -P . -M ${popMap}150 -n ${n_PARAM} -p 200  --catalog ./catalogQ2_125/catalog
# Move catalog to catalog dir.
mv catalog*gz ./catalogQ2_150/

#----------------- Rerun cstacks with another 25 samples ----------------------------
#Seventh round
time cstacks -P . -M ${popMap}175 -n ${n_PARAM} -p 200  --catalog ./catalogQ2_150/catalog
# Move catalog to catalog dir.
mv catalog*gz ./catalogQ2_175/

#----------------- Rerun cstacks with another 25 samples ----------------------------
#Last round
time cstacks -P . -M ${popMap}200 -n ${n_PARAM} -p 200  --catalog ./catalogQ2_175/catalog
# Move catalog to catalog dir.
mv catalog*gz ./catalogQ2_200/


echo "Finished the cstacks runs for 25, 50, 75, 100, 125 and 200 samples as input."


