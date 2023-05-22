#!/projects/mgdesaix@colostate.edu/mambaforge/envs/msprime/bin/python

import numpy as np
import sys

z_filename = sys.argv[1]
reference_filename = sys.argv[2]
outname = sys.argv[3]

z=np.loadtxt(z_filename)
ref=np.loadtxt(reference_filename)

ref_mean=np.mean(ref)
ref_std=np.std(ref)
new_z=(z - ref_mean)/ref_std

np.savetxt(outname, new_z, fmt="%.7f")
