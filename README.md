# kscale processing

Initial code to analyse results from K-Scale UM simulations. In particular, this includes comparison with satellite observations (CERES, will also add MIMIC). Using xarray as it is not fussy about exact data format. 

Files included:


Requires combined observation data in right format (scripts to be added)

Plans:
- Needs a bit more work to optimise performance and memory usage. In particular, will try to move to timestep-by-timestep processing wherever possible
- Generalise

Notes
- The simulations have some issues around metadata. In particular, it is hard to exactly align grids between different variables, and output times between variables and data sources. Sometimes, the initial or final dump are missing. To ensure the climatology is still reasonable, we first calculate a mean by hour-of-day, and then another mean over these hours to compare against e.g. CERES.
