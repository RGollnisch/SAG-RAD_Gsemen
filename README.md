Stacks runs for the SAG-RAD two continents project

The repo contains the scripts used to run the analyses for the SAG-RAD *Gonyostomum semen* project on 25 populations thoughout Europe and US.

In order to be able to evaluate the different settings of each of the steps in Stacks, the available pipeline denovo_map.pl was not used. 


**Before using any of the scripts below, please check that denovo_map.pl or ref_map.pl is not a better option.**

## Software requirements:
  Installed using Conda in PopGen environment:   
  
  **Kraken2** version 2.1.2 
  **Structure** version 2.3.4
  
  Installed separately either because the software was not available in Conda or had old version:   
  
  **Stacks** version 2.5.9
  For more information about running Stacks and the parameters used, please see the Stacks manual https://catchenlab.life.illinois.edu/stacks/manual/
  
## Project file structure
All analyses outputs are placed in separate sub directories under Analysis/ and scripts are expected to be run from the Analysis directory:
>   
├── Analysis   
│   .├── fastqc   
│   .├── kraken2   
│   .├── multiqc   
│   .├── ustacks   
│   .├── cstacks   
│   .├── sstacks   
│   .├── tsv2bam   
│   .├── populations   
│   .└── structure   



## Analyses performed

The following was done before this set of scripts were used:

0. Demultiplexing and cleaning of the fastq files was run on Uppmax.


After transferring to the local server the following analyses were performed:

## Contamination filtering
The cleaned and demultiplexed reads where filtered for contaminations of Human & Bacterial DNA using Kraken2 version 2.1.2. Any read pairs matching the databases were removed and the unclassified read where renamed to fit the required file naming of Stacks for PCR duplicate removal.

## ustacks 
Run ustacks.sh and change the values of:
	M_VALUE=4
	LT_m_VALUE=4
M_VALUE corresponds to the -M option in ustacks (Max distance between stacks)
LT_m_VALUE corresponds to the -m option (Min depth of coverage)
-N was kept at default values, which equals M+2 (Max distance allowed for adding second reads)

```{BASH}
nohup bash ustacks.sh
```

## cstacks
Run cstacks.sh 
After some initial tests it was decided to pick a representative set of samples from all populations according to the following criteria:
Select five samples from each population that belong to the upper half quartile of number of final stacks (min 50k loci) and have the highest mean coverage. 

A test of 50k loci was also performed.

This approach excludes some outliers that had a very high number of final stacks but low mean coverage, which are likely to be enriched for oversplit loci.

Note that the final number of individuals were lower than the 200 since some populations did not have enough samples that fulfilled the above mentioned criteria.

```{BASH}
nohup bash cstacks.sh
```

**Overview of RAD loci and Stacks catalog:**

## sstacks

Run sstacks.sh
Only the popmap (*50k* used for the final run) and catalog directory (catalog with 200 samples, see above) needs to be selected


```{BASH}
nohup bash sstacks.sh
```

## tsv2bam

At this step the renamed, filtered FASTQ reads are required in order to do PCR duplicate removal in gstacks. Note that the parameter -R points to a directory where all the FASTQ files should be located (not separated by lane anymore, hence the need to rename files). 

```{BASH}
nohup bash tsv2bam.sh
```

## gstacks
The paired end reads is checked for PCR duplicates with --rm-pcr-duplicates and the rest of the reverse reads is incorporated into an assembly of each contig.

```{BASH}
nohup bash gstacks.sh
```

## populations
For future simplification of testing multiple parameters it would be useful to have -r and -p to be set as bash parameters instead of being variables in the script itself.

Note that --write-single-snp -W  was not added to the populations since the loci remaining were below 1000.

E.g. --write-single-snp -W ./wl_1000

```{BASH}
nohup bash populations.sh
```


## Structure

In order to run Structure the following modifications of the Stacks structure output is required.

**1.** Remove the Stacks header line from the populations.stacks
Stacks adds a first line starting with #Stacks version...   
This needs to be removed before running Structure.

**2.** Check that the population names are integers: SM -> 1.
If the popmap provided in populations then the second column specifying the population must be changed into integers. E.g. Population SM was given the population id 1.

In order to effectively run Structure with many different parameters and also a minimum of five times (replicates with different random seeds), the parameters are modified on the command line and hence the same mainparams and extraparams can be used.

```{bash}

# Remove Stacks first line and change all USA population names to integers
 
# grep -v "^# Stacks" populations.structure | sed 's/SM /1  /'| sed 's/CC /2  /'|sed 's/CS /3  /'|sed 's/LL  /4  /'|sed 's/NO  /5  /'|sed 's/PQ  /6  /'|sed 's/PT  /7  /'|sed 's/RE  /8  /'|sed 's/VO  /9  /'|sed 's/WN  /10 /'> structure.EDT.txt

# Remove Stacks first line and change all EUR population names to integers
#Note that there should be tabs after each sample name & id. E.g. BR<tab> and 11<tab>.

#grep -v "^# Stacks" populations.structure | sed 's/BR /11 /'| sed 's/GF   /12 /'|sed 's/HE  /13 /'|sed 's/HU  /14 /'|sed 's/KO  /15 /'|sed 's/NT  /16 /'|sed 's/PB  /17 /'|sed 's/PE  /18 /'|sed 's/PI  /19 /'|sed 's/PS  /20 /'|sed 's/RM  /21 /'|sed 's/RO  /22 /'|sed 's/SE  /23 /'|sed 's/SI  /24 /'|sed 's/VP  /25 /'> structure.EDT.txt

# Run a loop to run one replicate for each K value
# Samples from USA has 10 populations, hence K 2-10 is run
# By changing the rand value in the -o parameter (rand1, rand2 ..) the output files will not overwrite previous replication.

# Check the number of variant sites retained (last value) and modify the -L value accordingly. 
#grep Kept populations.log

#Check the number of samples & populations and adjust -L for samples and K 2 to populations.
#grep Working populations.log

# variable r= Number of replications and variable k= the population sizes to be tested. Note that r*k jobs will be run in parallel.

# Note! The following script starts all structure runs in parallel. E.g. 10 repetitions * 14 K values= 140 processes!! Choose fewer k or r values to reduce this number if run on a smaller machine.

nohup bash structure.sh

```

**Analysis to determine the best value for K in Structure runs using the Evanno method:**
```
structure_evanno.R
```

## fineRADstructure
fineRADStructure was run following instrcutions on the fineRADstrucutre webiste: https://www.milan-malinsky.org/fineradstructure

## Fst figures and isolation by distance (IBD) analysis

Script to plot Fst matrix and isolation by distance (IBD) analysis:
```
fst_ibd_analysis.R
```

## Discriminant analysis of principal components (DAPC)

Discriminant analysis of principal components (DAPC) with environmental variables superimposed:
```
dapc.R
```

