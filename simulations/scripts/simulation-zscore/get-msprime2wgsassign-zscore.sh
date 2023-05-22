#!/bin/bash
#set a job name
#SBATCH --job-name=msprime2wgsassign-zscore
#SBATCH --output=./err-out/msprime2wgsassign-zscore.%A_%a.out
#SBATCH --error=./err-out/msprime2wgsassign-zscore.%A_%a.err
################
#SBATCH --time=24:00:00
#SBATCH --ntasks=20
#SBATCH --array=1-120
#################

source ~/.bashrc

for z in {1..5}
do

if [[ $SLURM_ARRAY_TASK_ID -eq 1 ]] && [[ $z -eq 1 ]]; then echo "sim" "rep" "samples" "read_depth" "fst_13" "fst_23" "ne_pop2" "acc_pop32" "acc_pop2" "z_ind_pop3" "pop" "runtime_ref_z" "runtime_z" >> /home/mgdesaix@colostate.edu/scratch/simulations/simulations-zscore/summary/simulation-zscore-summary.txt; fi

cd /home/mgdesaix@colostate.edu/scratch/simulations/simulations-zscore

migration=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $1}' migration-samples-readDepth-array-zscore.txt)
# ex. 0.001
samples=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $2}' migration-samples-readDepth-array-zscore.txt)
# ex. 10
read_depth=$(awk -v N=$SLURM_ARRAY_TASK_ID 'NR == N {print $3}' migration-samples-readDepth-array-zscore.txt)
# ex. 1

# ************************************************************
### Run simulation
# ************************************************************
conda activate msprime
full_vcf=$(./stepping-stone-pop-3.py ${migration} ${samples} ${read_depth} ${SLURM_ARRAY_TASK_ID} ${z})
# ex. ./msprime_vcfs/simulated_10_samples_0.01534-0.01433-0.05734_fst.gl_0.5_100_2.vcf

prefix=$(echo ${full_vcf} | cut -f3 -d/ | sed 's/.vcf//')
# ex. simulated_10_samples_0.01534_0.01433_0.05734_fst.gl_0.5_100_2.vcf
fst_13=$(echo ${prefix} | cut -f6 -d_)
fst_23=$(echo ${prefix} | cut -f5 -d_)

# ************************************************************
### Subset 100k and filter maf 100k
# ************************************************************
conda activate vcf
# remove underscore from vcf because of plink
sed -i 's/_//g' ${full_vcf}

# remove funky genotypes from high mutation rate (I think is the issue)
plink --vcf ${full_vcf} --recode --geno 0.0 --out ./msprime_vcfs/${prefix}

# LD filtering
plink --file ./msprime_vcfs/${prefix} --indep-pairwise 100 5 0.5 --allow-extra-chr --out ./msprime_vcfs/${prefix}.ld --maf 0.05
shuf -n 100000 ./msprime_vcfs/${prefix}.ld.prune.in > ./msprime_vcfs/${prefix}.filtered_100k
plink --file ./msprime_vcfs/${prefix} --extract ./msprime_vcfs/${prefix}.filtered_100k --out ./msprime_vcfs/${prefix}.filtered_100k --recode vcf
rm ./msprime_vcfs/${prefix}.ld.*

filtered=${prefix}.filtered_100k

rm ${full_vcf}

conda activate vcf

outname=${filtered}
maf=$(echo ${outname} | sed "s/$prefix\.//" | cut -f1 -d_)

# ************************************************************
### Simulate read depths
# ************************************************************
vcfgl=/projects/mgdesaix@colostate.edu/programs/vcfgl_ad/vcfgl
## Sample read depths from vcf to get GL
${vcfgl} -in ./msprime_vcfs/${outname}.vcf -out ./bcf/${outname} -err 0.01 -depth ${read_depth} -pos0 1

# ************************************************************
### Reference and unknown Beagle and allele depths
### Prep reference pops and "unknown" pop
# reference pops are the first 2 sets of individuals in the beagle
# "unknown" pop to test z-score is last individuals, number of samples
# ************************************************************
reference=${outname}.reference
unknown=${outname}.unknown
pop2=${outname}.pop2

bcftools query -l ./bcf/${outname}.bcf | head -n $((${samples}*2)) > ./bcf/${reference}.vcf.IDs
bcftools query -l ./bcf/${outname}.bcf | head -n $((${samples}*2 + 10)) | tail -n 10  > ./bcf/${pop2}.vcf.IDs
bcftools query -l ./bcf/${outname}.bcf | tail -n 10 > ./bcf/${unknown}.vcf.IDs


# get allele depths

bcftools query -f '[%AD{0}\t%AD{1}\t]\n' -S ./bcf/${reference}.vcf.IDs ./bcf/${outname}.bcf | gzip -c > ./counts/${reference}.counts.gz
bcftools query -f '[%AD{0}\t%AD{1}\t]\n' -S ./bcf/${pop2}.vcf.IDs ./bcf/${outname}.bcf | gzip -c > ./counts/${pop2}.counts.gz
bcftools query -f '[%AD{0}\t%AD{1}\t]\n' -S ./bcf/${unknown}.vcf.IDs ./bcf/${outname}.bcf | gzip -c > ./counts/${unknown}.counts.gz

## Scale GLs with awk code

bcftools query -f '%CHROM:%POS\t%REF\t%ALT{0}[\t%GL{0},%GL{1},%GL{2}]\n' -S ./bcf/${reference}.vcf.IDs ./bcf/${outname}.bcf | \
    sed 's/\.,\.,\./0,0,0/g' | \
    awk -f ./pl2gl.awk > ./beagles/${reference}.tmp

bcftools query -f '%CHROM:%POS\t%REF\t%ALT{0}[\t%GL{0},%GL{1},%GL{2}]\n' -S ./bcf/${pop2}.vcf.IDs ./bcf/${outname}.bcf | \
    sed 's/\.,\.,\./0,0,0/g' | \
    awk -f ./pl2gl.awk > ./beagles/${pop2}.tmp

bcftools query -f '%CHROM:%POS\t%REF\t%ALT{0}[\t%GL{0},%GL{1},%GL{2}]\n' -S ./bcf/${unknown}.vcf.IDs ./bcf/${outname}.bcf | \
    sed 's/\.,\.,\./0,0,0/g' | \
    awk -f ./pl2gl.awk > ./beagles/${unknown}.tmp

cd ./beagles

### Reference
n_col=$(awk '{print NF; exit}' ${reference}.tmp)
n_ind=$(echo $((n_col/3)))
snps=$(wc -l ${reference}.tmp | cut -f1 -d' ')

for ((i=1; i<=$n_ind; i++)); do if [[ $i -eq 1 ]]; then echo -ne "marker\tallele1\tallele2\t" >> ${reference}.header.txt; elif [[ $i -gt 1 ]] && [[ $i -lt $n_ind ]]; then ind=$(echo $((i-2))); echo -ne "Ind"${ind}"\tInd"${ind}"\tInd"${ind}"\t" >> ${reference}.header.txt;  else ind=$(echo $((i-2))); echo -ne "Ind"${ind}"\tInd"${ind}"\tInd"${ind}"\n" >> ${reference}.header.txt; fi; done

cat ${reference}.header.txt ${reference}.tmp | gzip -c > ${reference}.beagle.gz

rm ${reference}.header.txt
rm ${reference}.tmp

### Pop2 
n_col=$(awk '{print NF; exit}' ${pop2}.tmp)
n_ind=$(echo $((n_col/3)))

for ((i=1; i<=$n_ind; i++)); do if [[ $i -eq 1 ]]; then echo -ne "marker\tallele1\tallele2\t" >> ${pop2}.header.txt; elif [[ $i -gt 1 ]] && [[ $i -lt $n_ind ]]; then ind=$(echo $((i-2))); echo -ne "Ind"${ind}"\tInd"${ind}"\tInd"${ind}"\t" >> ${pop2}.header.txt; else ind=$(echo $((i-2))); echo -ne "Ind"${ind}"\tInd"${ind}"\tInd"${ind}"\n" >> ${pop2}.header.txt; fi; done

cat ${pop2}.header.txt ${pop2}.tmp | gzip -c > ${pop2}.beagle.gz

rm ${pop2}.header.txt
rm ${pop2}.tmp

# Unknown
n_col=$(awk '{print NF; exit}' ${unknown}.tmp)
n_ind=$(echo $((n_col/3)))

for ((i=1; i<=$n_ind; i++)); do if [[ $i -eq 1 ]]; then echo -ne "marker\tallele1\tallele2\t" >> ${unknown}.header.txt; elif [[ $i -gt 1 ]] && [[ $i -lt $n_ind ]]; then ind=$(echo $((i-2))); echo -ne "Ind"${ind}"\tInd"${ind}"\tInd"${ind}"\t" >> ${unknown}.header.txt;  else ind=$(echo $((i-2))); echo -ne "Ind"${ind}"\tInd"${ind}"\tInd"${ind}"\n" >> ${unknown}.header.txt; fi; done

cat ${unknown}.header.txt ${unknown}.tmp | gzip -c > ${unknown}.beagle.gz

rm ${unknown}.header.txt
rm ${unknown}.tmp

# ************************************************************
### WGSassign
# ************************************************************
conda activate WGSassign

cd ./WGSassign

#### Reference population stuff
for (( j=1; j<=$samples; j++ )); do echo $j pop1 >> ${reference}.IDs.txt; done
for (( j=1; j<=$samples; j++ )); do echo $j pop2 >> ${reference}.IDs.txt; done
sed -i 's/ /\t/g' ${reference}.IDs.txt

WGSassign --beagle ../${reference}.beagle.gz --pop_af_IDs ${reference}.IDs.txt --get_reference_af --ne_obs --out ${reference} --threads 20

ne_obs=${reference}.ne_obs.npy
mean_ne=$(./mean_ne_obs.py ${ne_obs})
# ne_pop1=$(echo ${mean_ne} | cut -f1 -d' ')
ne_pop2=$(echo ${mean_ne} | cut -f2 -d' ')

# Get reference z-score
start=`date +%s`
WGSassign --beagle ../${reference}.beagle.gz --pop_af_IDs ${reference}.IDs.txt --pop_names ${reference}.pop_names.txt --ind_ad_file ../../counts/${reference}.counts.gz --get_reference_z_score --out ${reference} --threads 20
end=`date +%s`
runtime_ref_z=$((end-start))

# Assign pop2 individuals
WGSassign --beagle ../${pop2}.beagle.gz --pop_af_file ${reference}.pop_af.npy --get_pop_like --out ${pop2} --threads 20
acc_pop2=$(awk '$2 > $1 {count++} END {print count/10}' ${pop2}.pop_like.txt)

# Pop2 z-scores
for (( j=1; j<=10; j++ )); do echo $j pop2 >> ${pop2}.IDs.txt; done
sed -i 's/ /\t/g' ${pop2}.IDs.txt
WGSassign --beagle ../${pop2}.beagle.gz --pop_af_IDs ${pop2}.IDs.txt --pop_af_file ${reference}.pop_af.npy --pop_names ${reference}.pop_names.txt --ind_ad_file ../../counts/${pop2}.counts.gz --get_assignment_z_score --out ${pop2} --threads 20

# Assign unknown individuals
WGSassign --beagle ../${unknown}.beagle.gz --pop_af_file ${reference}.pop_af.npy --get_pop_like --out ${unknown} --threads 20
acc_pop32=$(awk '$2 > $1 {count++} END {print count/10}' ${unknown}.pop_like.txt)

# Unknown z-scores
for (( j=1; j<=10; j++ )); do echo $j pop2 >> ${unknown}.IDs.txt; done
sed -i 's/ /\t/g' ${unknown}.IDs.txt
start=`date +%s`
WGSassign --beagle ../${unknown}.beagle.gz --pop_af_IDs ${unknown}.IDs.txt --pop_af_file ${reference}.pop_af.npy --pop_names ${reference}.pop_names.txt --ind_ad_file ../../counts/${unknown}.counts.gz --get_assignment_z_score --out ${unknown} --threads 20
end=`date +%s`
runtime_z=$((end-start))

tail -n ${samples} ${reference}.reference_z_ind.txt > ${reference}.reference_z_ind.pop2.txt


./standardize_z.py ${unknown}.z_ind.txt ${reference}.reference_z_ind.pop2.txt ${unknown}.z_ind.standardized.txt
./standardize_z.py ${pop2}.z_ind.txt ${reference}.reference_z_ind.pop2.txt ${pop2}.z_ind.standardized.txt

for z_ind_pop3 in `cat ${unknown}.z_ind.standardized.txt`; do echo ${SLURM_ARRAY_TASK_ID} ${z} ${samples} ${read_depth} ${fst_13} ${fst_23} ${ne_pop2} ${acc_pop32} ${acc_pop2} ${z_ind_pop3} "pop3" ${runtime_ref_z} ${runtime_z} >> /home/mgdesaix@colostate.edu/scratch/simulations/simulations-zscore/summary/simulation-zscore-summary.txt; done
# header = "sim" "rep" "samples" "read_depth" "fst_13" "fst_23" "ne_pop2" "acc_pop32" "acc_pop2" "z_ind_pop3" "pop" "runtime_ref_z" "runtime_z"

for z_ind_pop3 in `cat ${pop2}.z_ind.standardized.txt`; do echo ${SLURM_ARRAY_TASK_ID} ${z} ${samples} ${read_depth} ${fst_13} ${fst_23} ${ne_pop2} ${acc_pop32} ${acc_pop2} ${z_ind_pop3} "pop2" ${runtime_ref_z} ${runtime_z} >> /home/mgdesaix@colostate.edu/scratch/simulations/simulations-zscore/summary/simulation-zscore-summary.txt; done

rm ${outname}.*
rm ../${outname}.*
rm ../../bcf/${outname}.*
rm ../../msprime_vcfs/${outname}.*
rm ../../msprime_vcfs/${prefix}.*
rm ../../counts/${outname}.*

done

