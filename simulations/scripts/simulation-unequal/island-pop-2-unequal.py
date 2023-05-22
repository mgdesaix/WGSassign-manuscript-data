#!/projects/mgdesaix@colostate.edu/mambaforge/envs/msprime/bin/python

import numpy as np
import msprime
import sys

migration = float(sys.argv[1])
samples1 = int(sys.argv[2])
samples2 = int(sys.argv[3])
read_depth = float(sys.argv[4])
slurm_id = int(sys.argv[5])
iteration = int(sys.argv[6])

demography = msprime.Demography.island_model([1000, 1000], migration_rate = migration)
ts = msprime.sim_ancestry({0:samples1, 1:samples2},
                          demography = demography,
                          recombination_rate=1e-8,
                          sequence_length=1e8)

mts = msprime.sim_mutations(ts, rate = 7e-7)

## Calculate Fst
A = mts.samples()[:samples1*2]
B = mts.samples()[samples1*2:]
d = mts.Fst([A,B])
fst = round(d, 5)

## Save VCF

outname = "./msprime_vcfs/simulated_" + str(samples1) + "_" + str(samples2) + "_samples_" + str(fst) + "_fst.gl_" + str(read_depth) + "_" + str(slurm_id) + "_" + str(iteration) + ".vcf"

with open(outname, "w") as vcf_file:
    mts.write_vcf(vcf_file)

print(outname)
