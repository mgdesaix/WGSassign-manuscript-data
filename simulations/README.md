Genetic simulations for WGSassign
================
Matt G. DeSaix


In the WGSassign manuscript, we implemented 3 genetic simulations to assess the behavior of the WGSassign population assignment tests:

1.  Two-population island model with equal sample sizes
2.  Two-population island model with unequal sample sizes
3.  Three-population stepping-stone model for assessing z-score metric

All of these simulations followed similar general workflows:

* Perform genetic simulations in [msprime](https://tskit.dev/msprime/docs/stable/intro.html)
* Filter variants with [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* Simulate genotype likelihoods from given read depths in [vcfgl](https://github.com/isinaltinkaya/vcfgl)
* Convert VCF files to Beagle files (custom script: pl2gl.awk)
* Run WGSassign

The intermediary output from all of these steps is too large to provide the raw data. However, the summary data from these three simulations is provided in the `data/` directory. All scripts are provided in the `scripts/` directory, organized by simulation except for scripts that remained the same for all analyses. The `.sh` files performed the entire workflow from simulating genetic data to running WGSassign. The `.Rmd` file was used to summarize all data and make figures for the simulated data.

