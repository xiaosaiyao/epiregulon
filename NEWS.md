# epiregulon 2.0.1
Fixed bugs that caused `addWeights` to fail when `aggregateCells = TRUE`.
When `aggregateCells` is set to `TRUE`, cells are aggregated by calculating the mean of the original cell features rather than their sums.
Fixed a warning about the class of the `cellNum` argument passed to `calculateP2G`, which was raised incorrectly.


# epiregulon 2.0.0
* New function in the workflow, `optimizeMatacellNumber`, which estimates the optimal value of the `cellNum` parameter
passed to `calculateP2G`.
* Using `scrapper::clusterKmeans` to implement deterministic algorithm of determination k-mean cluster.
It is used to create metacells in the `calculateP2G` function.
* Added to the output of `calculateP2G` function empirical p-values and FDR. The null distribution
is calculated based on random links of peaks to the genes from other chromosomes.
* Default ChIP-seq output from `tfBinding` is changed to version 2 which imposes more stringent cutoffs.
We have increased the number of unique reads to 20M (previously 10M in version 1) and the number of peaks passing 
FDR < 1e-5 to 1000 peaks (previously 100 peaks in version 1). For ENCODE data, we now remove any samples with 
Audit.NOT_COMPLIANT or Audit.ERROR. Flags.
* The argument `method` to `calculateActivity` has been deprecated.
* Bug corrected in `addLogFC`. p-values were mistakenly transformed from log10 when it should be transformed from natural log 
when `logFC_ref` or/and `logFC_condition` was specified 

# epiregulon 1.5.1
* checking for the duplicated gene names in the input gene expression SingleCellExperiment
* validation of the regulon object passed to pruneRegulon, addWeights, and calculateActivity 

# epiregulon 1.0.1
* the function `addMotifScore` returns correct values instead of NA.

# epiregulon 0.99.20
* the function `aggregateAcrossCellsFast` has been added.

# epiregulon 0.9.12
* in `calculateP2G` staring position of the gene has been corrected.

# epiregulon 0.99.11
* new function `aggregateAcrossCells` has been added, which is implemented in c++ and makes 
run time for `calculateP2G` much shorter. 

# epiregulon 0.99.5

* new function `addLogFC` which adds log fold changes of gene expression to regulons and significance statistics for differential gene expression
* `cellNum` argument to `calculateP2G` set by default to 100 (previously it was 200).
* in `pruneRegulon` the check to the uniqueness of gene names has been added.

# epiregulon 0.99.0
Version number downgraded to 0.99.0 to meet Bioconductor requirements.

# epiregulon 0.98.0

commit 12fabd3278b76711f4b2cbdcaca6794a918ffd83

* transcription factor motif databeses used by `addMotifScore` function (`human_pwms_v1`, and `mouse_pwms_v1`) were replaced wih the new ones (`human_pwms_v2`, `mouse_pwms_v2`) which are up-to-date.

BUG FIX

* `clusters` vector provided to the `addWeights` and `pruneRegulon` functions is checked for the presence of NA values. If they are present the function stops preventing segmentation fault and session abortion.

