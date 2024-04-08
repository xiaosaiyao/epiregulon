# epiregulon 0.99.11

# epiregulon 0.99.5

* new function `addLogFC` which adds log fold changes of gene expression to regulons and significance statistics for differential gene expression
* `cellNum` argument to `calculateP2G` set by default to 100 (previously it was 200).
* in `pruneRegulon` the check to the uniqueness of gene names has been added.

# epiregulon 0.99.0
Version number downgraded to 0.99.0 to meet Bioconductor requirements.

# epiregulon 1.0.38

* transcription factor motif databeses used by `addMotifScore` function (`human_pwms_v1`, and `mouse_pwms_v1`) were replaced wih the new ones (`human_pwms_v2`, `mouse_pwms_v2`) which are up-to-date.

BUG FIX

* `clusters` vector provided to the `addWeights` and `pruneRegulon` functions is checked for the presence of NA values. If they are present the function stops preventing segmentation fault and session abortion.

