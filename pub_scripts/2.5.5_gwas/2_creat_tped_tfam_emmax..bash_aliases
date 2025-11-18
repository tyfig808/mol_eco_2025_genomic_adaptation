#!/usr/bin/env bash
#SBATCH --account=schiestl.systbot.uzh
#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --time=24:00:00
#SBATCH --array=0-7%8

module load anaconda3
eval "$(conda shell.bash hook)"
source activate vcftools

wd=~/mol_eco_2024/2.5.5_gwas

f=("$wd"*"G1.numchr.bed")
i=$SLURM_ARRAY_TASK_ID

number1=$(echo ${f[i]} | awk -F / '{ print $NF}')
number2=$(echo $number1 | awk -F . '{ print $1 "." $2}')

echo "$number2"
### create tped and tfam for gwas with emmax
plink2 --bfile $wd"$number2" --allow-extra-chr --recode 12 transpose --maf 0.05 --remove remove_file --output-missing-genotype 0 --out $wd"$number2"