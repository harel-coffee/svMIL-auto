#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=4G
#$ -l h_rt=00:10:00
#$ -e eQTLCorrelation_err
#$ -o eQTLCorrelation_out

#Run 1 normal run without permutations

#Then run the permutation runs 1000 times, create jobs for that. 


uuid=$(uuidgen)

permutationsYN="False"

qsub determineCorrelatedEQTLPairs.sh "$permutationsYN" "$uuid" 0 #dummy value of 0 because we do not use permutation numbers here

permutationsYN="True"

qsub -t 1-100:1 -tc 100 determineCorrelatedEQTLPairs.sh "$permutationsYN" "$uuid" #use job array and inside this script SGE TASK ID to get the run number
