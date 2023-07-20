# CNN_Canarias

scripts from Rincon-Barrado et al. "Micro and macroevolution in the Rand Flora pattern: new advances in analytical tools integrating phylogenomics and phylogeography", to simulate and compare genetic datasets of two species from the Canary islands under different demographic scenarios. The selected model for each species was then used to perform parameter estimation. Both model comparison and parameter estimation were achieved by training a CNN, from which prediction were then used as Summary Statistics (SuSt) for ABC.

simulate_Euphorbia.py - python script to simulate segregating sites and save them as NumPy arrays for *Euphorbia balsamifera*. Simulations for the other species use a similar script, modified to reflect differences in the sample sizes and priors.

TrainCalibrate_ModSel.py - python script for training and calibrating a CNN for model selection, perform cross-validation and predict the most likely model from empirical data in the species *E. balsamifera* (script for *Kleinia* is similar).

TrainParameterEst.py - python script for training a CNN to peform parameter estimation in the species *E. balsamifera* (script for *Kleinia* is similar).

ABC_Canarias.R - R script for ABC cross-validation and empirical data prediction for model selection and parameter estimation in the species *Euphorbia balsamifera* (script for *Kleinia* is similar).
