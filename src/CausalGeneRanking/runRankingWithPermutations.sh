#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=4G
#$ -l h_rt=00:10:00
#$ -e ranking_permutations_err
#$ -o ranking_permutations_out

#Run 1 normal run without permutations

#Then run the permutation runs 1000 times, create jobs for that. 

causalGeneFile="$1"

uuid=$(uuidgen)

permutationsYN="False"

mode="$2"

qsub runRanking.sh "$causalGeneFile" "$uuid" "$permutationsYN" "$mode"

permutationsYN="True"

qsub -t 1-1000:1 -tc 100 runRanking.sh "$causalGeneFile" "$uuid" "$permutationsYN" "$mode" #use job array and inside this script SGE TASK ID to get the run number
