Genetic simulations for WGSassign
================
Matt G. DeSaix


In the WGSassign manuscript, we implemented 3 genetic simulations to assess the behavior of the WGSassign population assignment tests:

1.  Two-population island model with equal sample sizes
2.  Two-population island model with unequal sample sizes
3.  Three-population stepping-stone model for assessing z-score metric

All of these simulations followed similar general workflows:

* Perform genetic simulations in [msprime](https://tskit.dev/msprime/docs/stable/intro.html)
* Simulate genotype likelihoods from given read depths in [vcfgl](https://github.com/isinaltinkaya/vcfgl)
* Filter variants with [bcftools](https://samtools.github.io/bcftools/bcftools.html)
* Convert VCF files to Beagle files (custom script)
* Run WGSassign

