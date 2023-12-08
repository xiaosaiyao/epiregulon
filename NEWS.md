# epiregulon 1.0.38

* transcription factor motif databeses used by `addMotifScore` function (`human_pwms_v1`, and `mouse_pwms_v1`) were replaced wih the new ones (`human_pwms_v2`, `mouse_pwms_v2`) which are up-to-date.

BUG FIX

* `clusters` vector provided to the `addWeights` and `pruneRegulon` functions is checked for the presence of NA values. If they are present the function stops preventing segmentation fault and session abortion.
