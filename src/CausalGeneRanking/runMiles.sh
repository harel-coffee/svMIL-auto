#!/bin/bash
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -m as
#$ -M m.m.nieboer@umcutrecht.nl
#$ -l h_vmem=50G
#$ -l h_rt=24:00:00
#$ -e miles_lasso_SG_pr_err
#$ -o miles_lasso_SG_pr_out

python milesAllBags.py