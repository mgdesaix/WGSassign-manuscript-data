Yellow warbler and WGSassign
================
Marina D. Rodriguez

This repository contains the code used to empirically test population assignment in WGSassign on Yellow Warbler data.

The beagle file and ID file needed as input can be accessed through Dryad: 

To measure the assignment accuracy of WGSassign, we used leave-one-out cross validation (the --loo specification in WGSassign) using the input beagle file (yewa.known.ind105.ds_2x.beagle.gz) and our ID file (yewa.known.ind105.reference.IDs.txt). The ID file is a tab-delimited file with 2 columns, the first being the sample ID, and the second being the known reference population. The sample order in the ID file should match that of the input beagle file. Leave-one-out cross validation was run using the script [WGSassign.loo.sh](https://github.com/mgdesaix/WGSassign-manuscript-data/blob/main/yellow-warbler/WGSassign.loo.sh).

We then checked the assignment accuracy using the output from leave-one-out cross validation in R. Using the provided R-script ([WGS_assignment.check_yewa.Rmd](https://github.com/mgdesaix/WGSassign-manuscript-data/blob/main/yellow-warbler/WGSA_assignment.check_yewa.Rmd)), we first read in our ID file and our output from WGSassign (yewa.wgs.assignment_LOO.txt). The provided script was used to then summarize and provide population-specific assignment accuracy.
