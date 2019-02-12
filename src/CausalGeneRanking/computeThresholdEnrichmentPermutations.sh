#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=4G
#$ -l h_rt=00:30:00
#$ -e threshold_p_err
#$ -o threshold_p_out

#Create jobs to compute threshold enrichment on every permutation file

fileNum="$SGE_TASK_ID"

python computeThresholdEnrichment.py "$fileNum" "../../data/Genes/Census_allTue\ Apr\ 10\ 14_56_44\ 2018.tsv" "../../data/Genes/genesWithSNVs.txt" "../../data/DEG/BRCA.txt" "RankedGenes/69d48943-7dca-4df1-b6a6-13f6f78e83bc/BRCA/"
