#!/bin/bash
#set a job name
#SBATCH --job-name=msprime2wgsassign-maf
#SBATCH --output=./err-out/msprime2wgsassign-maf.%A_%a.out
#SBATCH --error=./err-out/msprime2wgsassign-maf.%A_%a.err
################
#SBATCH --time=24:00:00
#SBATCH --ntasks=20
#SBATCH --array=1-480
#################

source ~/.bashrc

for z in {1..5}
do

if [[ $SLURM_ARRAY_TASK_ID -eq 1 ]] && [[ $z -eq 1 ]]; then echo "sim" "rep" "samples" "read_depth" "fst" "ne_pop1" "acc_pop1" "ne_pop2" "acc_pop2" "runtime" "runtime_full" "maf" "snps" >> /home/mgdesaix@colostate.edu/scratch/simulations/summary/simulation-maf-summary.txt; fi

cd /home/mgdesaix@colostate.edu/scratch/simulations

migration=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' migration-samples-readDepth-array-full.txt)
# ex. 0.001
samples=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' migration-samples-readDepth-array-full.txt)
# ex. 10
read_depth=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $3}' migration-samples-readDepth-array-full.txt)
# ex. 1

# ************************************************************
### Run simulation
# ************************************************************
conda activate msprime
full_vcf=$(./island-pop-2-equal.py ${migration} ${samples} ${read_depth} ${SLURM_ARRAY_TASK_ID} ${z})
# ex. ./msprime_vcfs/simulated_10_samples_0.01534_fst.gl_0.5_100_2.vcf

prefix=$(echo ${full_vcf} | cut -f3 -d/ | sed 's/.vcf//')
# ex. simulated_10_samples_0.01534_fst.gl_0.5_100_2
fst=$(echo ${prefix} | cut -f4 -d_)

# ************************************************************
### Subset 100k and filter maf 100k
# ************************************************************
conda activate vcf

# remove underscore from vcf because of plink
sed -i 's/_//g' ${full_vcf}
# remove funky genotypes from high mutation rate (I think is the issue)
plink --vcf ${full_vcf} --recode --geno 0.0 --out ./msprime_vcfs/${prefix}

# get raw IDs from plink for shuffling
cut -f2 ./msprime_vcfs/${prefix}.map | shuf -n 100000 > ./msprime_vcfs/${prefix}.raw_100k
plink --file ./msprime_vcfs/${prefix} --extract ./msprime_vcfs/${prefix}.raw_100k --out ./msprime_vcfs/${prefix}.raw_100k --recode vcf

# LD filtering
plink --file ./msprime_vcfs/${prefix} --indep-pairwise 100 5 0.5 --allow-extra-chr --out ./msprime_vcfs/${prefix}.ld --maf 0.05
shuf -n 100000 ./msprime_vcfs/${prefix}.ld.prune.in > ./msprime_vcfs/${prefix}.filtered_100k
plink --file ./msprime_vcfs/${prefix} --extract ./msprime_vcfs/${prefix}.filtered_100k --out ./msprime_vcfs/${prefix}.filtered_100k --recode vcf
rm ./msprime_vcfs/${prefix}.ld.*

raw=${prefix}.raw_100k
filtered=${prefix}.filtered_100k

rm ${full_vcf}

# go thru following steps for both the raw 100k and filtered 100k vcfs
for y in ${raw} ${filtered}
do
conda activate vcf
cd /home/mgdesaix@colostate.edu/scratch/simulations

outname=${y}
maf=$(echo ${outname} | sed "s/$prefix\.//" | cut -f1 -d_)

# ************************************************************
### Simulate read depths
# ************************************************************
vcfgl=/projects/mgdesaix@colostate.edu/programs/vcfgl/vcfgl
## Sample read depths from vcf to get GL
${vcfgl} -in ./msprime_vcfs/${outname}.vcf -out ./bcf/${outname} -err 0.01 -depth ${read_depth} -pos0 1

# ************************************************************
### Convert to Beagle
# ************************************************************
## Scale GLs with awk code

bcftools query -f '%CHROM:%POS\t%REF\t%ALT{0}[\t%GL{0},%GL{1},%GL{2}]\n' ./bcf/${outname}.bcf | \
    sed 's/\.,\.,\./0,0,0/g' | \
    awk -f ./pl2gl.awk > ./beagles/${outname}.tmp

cd ./beagles

n_col=$(awk '{print NF; exit}' ${outname}.tmp)
n_ind=$(echo $((n_col/3)))
snps=$(wc -l ${outname}.tmp | cut -f1 -d' ')

for ((i=1; i<=$n_ind; i++)); do if [[ $i -eq 1 ]]; then echo -ne "marker\tallele1\tallele2\t" >> ${outname}.header.txt; elif [[ $i -gt 1 ]] && [[ $i -lt $n_ind ]]; then ind=$(echo $((i-2))); echo -ne "Ind"${ind}"\tInd"${ind}"\tInd"${ind}"\t" >> ${outname}.header.txt;  else ind=$(echo $((i-2))); echo -ne "Ind"${ind}"\tInd"${ind}"\tInd"${ind}"\n" >> ${outname}.header.txt; fi; done

cat ${outname}.header.txt ${outname}.tmp | gzip -c > ${outname}.beagle.gz

rm ${outname}.header.txt
rm ${outname}.tmp

# ************************************************************
### WGSassign
# ************************************************************
conda activate WGSassign

cd ./WGSassign

for i in pop1 pop2; do for (( j=1; j<=$samples; j++ )); do echo $j $i >> ${outname}.IDs.txt; done; done
sed -i 's/ /\t/g' ${outname}.IDs.txt

start=`date +%s`
WGSassign --beagle ../${outname}.beagle.gz --pop_af_IDs ${outname}.IDs.txt --get_reference_af \
    --ne_obs --out ${outname} --threads 20
end=`date +%s`

runtime=$((end-start))

start=`date +%s`
WGSassign --beagle ../${outname}.beagle.gz --pop_af_IDs ${outname}.IDs.txt --get_reference_af --loo \
    --out ${outname} --threads 20
end=`date +%s`

runtime_full=$((end-start))

acc_pop1=$(awk -v N=$samples 'NR <= N && $1 > $2 {count++} END {print count/N}' ${outname}.pop_like_LOO.txt)
acc_pop2=$(awk -v N=$samples 'NR > N && $2 > $1 {count++} END {print count/N}' ${outname}.pop_like_LOO.txt)

ne_obs=${outname}.ne_obs.npy

mean_ne=$(./mean_ne_obs.py ${ne_obs})
ne_pop1=$(echo ${mean_ne} | cut -f1 -d' ')
ne_pop2=$(echo ${mean_ne} | cut -f2 -d' ')

echo ${SLURM_ARRAY_TASK_ID} ${z} ${samples} ${read_depth} ${fst} ${ne_pop1} ${acc_pop1} ${ne_pop2} ${acc_pop2} ${runtime} ${runtime_full} ${maf} ${snps} >> /home/mgdesaix@colostate.edu/scratch/simulations/summary/simulation-maf-summary.txt

rm ${outname}.*
rm ../${outname}.*
rm ../../bcf/${outname}.*
rm ../../msprime_vcfs/${outname}.*

done
rm /home/mgdesaix@colostate.edu/scratch/simulations/msprime_vcfs/${prefix}.*
done

