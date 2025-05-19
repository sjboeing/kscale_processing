# kscale processing

Initial code to evaluate UM K-Scale simulations. In particular, this includes comparison with satellite observations (CERES, will also add MIMIC). 

The scripts use xarray as this is not fussy about exact data format. 

Files included:
- dummy_config.yaml: example of a path/directory configuration with dummy filenames.
- dyamond_radiation_plots.py: plotting scripts and data types corresponding to plots.
- dyamond_setups.py: currently includes both the data types that describe the different dyamond runs and areas for detailed analysis, as well as the details of the current runs and the basins we analyse (which could be separated or brought into yaml files).
- radiation_process.py: currently includes all the processing scripts, could be separated out further.

Requires combined observation data in right format (scripts to be added).

Plans:
- Needs a bit more work to optimise performance and memory usage. In particular, it may be worth moving to timestep-by-timestep processing wherever possible
- Generalise

Notes
- The main challenge for this exercise is getting the metadata to be consistent. In particular, it is hard to exactly align grids between different variables and models, and to get consistent output times between model and observations. Sometimes, the initial or final dump are missing from the model. To ensure the climatology is still reasonable, we first calculate a mean by hour-of-day, and then another mean over these hours to compare against e.g. CERES.
