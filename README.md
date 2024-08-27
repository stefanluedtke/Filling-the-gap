# 🐟 Filling the Gap: Improving the spatio-temporal coverage of small pelagic fish surveys through modelling approaches

Scripts and datasets to reproduce the results of our manuscript on data imputation of fisheries-independent survey data. 
Specifically, we focus on acoustic surveys in the Baltic Sea aiming at estimating abundance and biomass of the pelagic species sprat and herring.
They are conducted biannually in May (Baltic Acoustic Spring Survey (BASS) resulting in an index for sprat) and October (Baltic International Acoustic survey 
(BIAS) resulting in an index for each herring and sprat. In this repository, we evaluate three  
modelling approaches: Linear mixed-effects models (LMMs), Gradient Boosted Trees (XGBoost) and Generalized 
Additive Models (GAMs) alongside the baseline approach. 

## Running the scripts

The script example.R shows how to use the code to load data, plot it and run imputation models. The analyses in the paper can be reproduced with the scripts in the folder eval. 
