#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=16000
#SBATCH --time=48:00:00
#SBATCH --array=2-9%8


# set wd and get files
wd=~/mol_eco_2024/2.5.3_sel_scans
files=("$wd"*"_g1_g10_pruned.sync")

# for bash array
i=$SLURM_ARRAY_TASK_ID

# read over ne file 
ne_file=$wd"avg_ne.csv"

# read line depending on slurm
line=$(awk 'NR=='$i'{ print; exit }' $ne_file)

# get treatment
treat=$(echo $line | awk -F , '{ print $1}')
echo "This is the vcf file being processed $treat ..."

# get ne
ne=$(echo $line | awk -F , '{ print $4}')
echo "This is the effective pop size for treatment $treat $ne ..."

#name sync file and out file
s_file='treatment_'$treat'_g1_g10_pruned.sync'
o_file='clear_'$treat'.out' 

# run clear 
python2 ~/CLEAR/CLEAR.py --sync $wd$s_file --N $ne --out $wd$o_file --plot

