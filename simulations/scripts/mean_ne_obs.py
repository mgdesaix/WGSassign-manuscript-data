#!/projects/mgdesaix@colostate.edu/mambaforge/envs/msprime/bin/python

import numpy as np
import sys

filename = sys.argv[1]

ne_obs=np.load(filename)
ne_means = np.mean(ne_obs, 0)

print(ne_means[0], ne_means[1])
