# MetCoOp-Ensemble-Prediction-System

This repository is used to download forecasts of weather variables from the MetCoOp Ensemble Prediction System (MEPS).

*download_meps.py: Downloads selected NWP variables for a chosen time period and lead times into the future.
*download_meps.py: Downloads wind speed at 10 meters and at a chosen hub height through the Hypsometric equation for a chosen time period and lead times into the future. In this way, the variability (e.g., diurnal) that different between the heights are taken into account.

For information about how this could be retrieved, see the Norwegian Meteorological Institute:
https://github.com/metno/NWPdocs/wiki

The datasets are found here:
https://thredds.met.no/thredds/metno.html

For information about the model parametrization and physics, see:

[1] Bengtsson, L. et al. (2017), "The HARMONIE–AROME Model Configuration in the ALADIN–HIRLAM NWP System", DOI: https://doi.org/10.1175/MWR-D-16-0417.1 
