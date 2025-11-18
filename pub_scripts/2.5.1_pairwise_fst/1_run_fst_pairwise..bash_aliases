module load anaconda3
eval "$(conda shell.bash hook)"
source activate vcftools

VCF=~/final_filter.vcf.gz
wd=~/fst/
o=~/fst/rep_sep/
f=("$wd"*"_rep_pop")
i=$SLURM_ARRAY_TASK_ID
for ((j = i + 1; j < ${#f[@]}; j++)); do
    number1=$(echo ${f[i]} | awk -F / '{ print $NF}')
    number2=$(echo ${f[j]} | awk -F / '{ print $NF}')
    out="$o$number1$number2"
    echo "$number1 - $number2"
    vcftools --gzvcf "$VCF" --weir-fst-pop ${f[i]} --weir-fst-pop ${f[j]} --fst-window-size 5000 --out "$out"
done