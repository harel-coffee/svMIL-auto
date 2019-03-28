#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=8G
#$ -l h_rt=02:00:00
#$ -e simpleMil_e
#$ -o simpleMil_o

#Create jobs to compute threshold enrichment on every permutation file

fileNum="$SGE_TASK_ID"

python runSimpleMIL.py geneSVPairs_somaticBRCA.txt geneSVPairs_germline.txt bagMIL/ "$fileNum"
