#epiregulon 0.99.20
* the function aggregateAcrossCellsFast has been added.

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

