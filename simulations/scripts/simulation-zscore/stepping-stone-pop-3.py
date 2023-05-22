#!/projects/mgdesaix@colostate.edu/mambaforge/envs/msprime/bin/python

import numpy as np
import msprime
import sys

migration = float(sys.argv[1])
samples1 = int(sys.argv[2])
read_depth = float(sys.argv[3])
slurm_id = int(sys.argv[4])
iteration = int(sys.argv[5])

demography = msprime.Demography.stepping_stone_model([1000, 1000, 1000], migration_rate = migration, boundaries = True)
ts = msprime.sim_ancestry({0:samples1, 1:samples1+10, 2:samples1+10},
                          demography = demography,
                          recombination_rate=1e-8,
                          sequence_length=1e8)

mts = msprime.sim_mutations(ts, rate = 7e-7)

## Calculate Fst
A = mts.samples()[:samples1*2]
B = mts.samples()[samples1*2:samples1*4]
C = mts.samples()[samples1*4+20:samples1*6+20]
d_AB = round(mts.Fst([A,B]), 5)
d_BC = round(mts.Fst([C,B]), 5)
d_AC = round(mts.Fst([A,C]), 5)

## Save VCF
# ex. ./msprime_vcfs/simulated_10_samples_0.01534_0.01433_0.05734_fst.gl_0.5_100_2.vcf
outname = "./msprime_vcfs/simulated_" + str(samples1) + "_samples_" + str(d_AB) + "_" + str(d_BC) + "_" + str(d_AC) + "_fst.gl_" + str(read_depth) + "_" + str(slurm_id) + "_" + str(iteration) + ".vcf"

with open(outname, "w") as vcf_file:
    mts.write_vcf(vcf_file)

print(outname)
