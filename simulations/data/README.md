Simulation data README
================

The data for the 3 simulations are provided in this directory. The column names for each file are as follows:

`simulation-equal-summary.txt`
There is one row of data summarizing each unique combination of `sim` and `rep`.

sim = simulation number unique to a given migration rate, read depth, and number of samples per population
rep = iteration of simulation workflow using the same parameters for each simulation number
samples = number of individuals sampled from each population
read_depth = mean read depth used to simulate genotype likelihoods
fst = Fst calculated between the two populations
ne_pop1 = effective sample size of population 1
acc_pop1 = leave-one-out assignment accuracy of individuals from population 1 (proportion of individuals correctly assigned)
ne_pop2 = effective sample size of population 2
acc_pop2 = leave-one-out assignment accuracy of individuals from population 2 (proportion of individuals correctly assigned)
runtime = time in seconds of allele frequency, fisher information, and effective sample size calculations
runtime_full = time in seconds of the above runtime + leave-one-out assignment
maf = raw (no variant filtering) or maf (variant filtering my minor allele frequency)
snps = number of SNPs as input into WGSassign

`simulation-equal-summary.txt`
There is one row of data summarizing each unique combination of `sim` and `rep`.
Column names are the same as above, except for:

samples_pop1 = number of individuals sampled from population 1
samples_pop2 = number of individuals sampled from population 2

`simulation-zscore-summary.txt`
There is one row of data for each individual.
Column names are the same as above, except for:

fst_13 = Fst calculated between population 1 and 3
fst_23 = Fst calculated between population 2 and 3
ne_pop2 = Effective sample size for population 2
acc_pop32 = Proportion of individuals from population 3 assigned to population 2
acc_pop2 = Proportion of individuals from population 2 assigned to population 2
z_ind = z-score for the individual
pop = Population of origin of the individual (population 2 or 3)
runtime_ref_z = Time in seconds taken for WGSassign to calculate the reference population z-scores
runtime_z = Time in seconds taken for WGSassign to calculate the z-scores for all individuals being assigned
