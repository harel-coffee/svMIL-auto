#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=16G
#$ -l h_rt=00:30:00
#$ -e eQTLCorrelation_normal_err
#$ -o eQTLCorrelation_normal_out

#Run 1 normal run without permutations

python determineEQTLCorrelatedPairs.py "$1" "$2" "$SGE_TASK_ID"
