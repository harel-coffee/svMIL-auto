#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=8G
#$ -l h_rt=05:00:00
#$ -e opt_err
#$ -o opt_out

python main.py "../Output/RankedGenes/milEnhProEqtlSe/BRCA/bags.pkl" "../pValues_allGenes_smallTads.txt"

