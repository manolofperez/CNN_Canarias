# CNN_Canarias

scripts from Rincon-Barrado et al. "Phylogenomic reconstruction and deep learning reveal an ancient Rand Flora origin and a recent admixture history in the Canarian multi-island Endemic Kleinia neriifolia", to simulate and compare genetic datasets of Kleinia neriifolia from the Canary islands under different demographic scenarios.

simulate_Euphorbia.py - python script to simulate segregating sites and save them as NumPy arrays for *Euphorbia balsamifera*. Simulations for the other species use a similar script, modified to reflect differences in the sample sizes and priors.

TrainCalibrate_ModSel.py - python script for training and calibrating a CNN for model selection, perform cross-validation and predict the most likely model from empirical data in the species *E. balsamifera* (script for *Kleinia* is similar).

TrainParameterEst.py - python script for training a CNN to peform parameter estimation in the species *E. balsamifera* (script for *Kleinia* is similar).

ABC_Canarias.R - R script for ABC cross-validation and empirical data prediction for model selection and parameter estimation in the species *Euphorbia balsamifera* (script for *Kleinia* is similar).
